#!/usr/bin/python
#==============================================================================
# calc-kpseudo-quant.py
#
# Shawn Driscoll
# 20170821
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Building a script that can combine kallisto pseudo output with the 
# isoform balance created in kallisto
#==============================================================================

import sys
import argparse
import math
import re
import os
import traceback
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	samples = []
	kpseudo = {}
	dannot = None
	dgid2tid = None
	sample_abundance = []

	# find the kalliso pseudo output. the 'cells' file in that output will have the sample names
	# which we can use to find all of the abundance.tsv files
	kp = KPseudo("kpseudo")
	if not kp.check_path():
		message("Error: kallisto pseudo output seems to be missing some files")
		return 1
	else:
		message("Kallisto pseudo output appears to be intact")
		# load it up
		kp.load_data()

	# load the fai
	lfai, tid2faiIndex = load_fai(args.fai)
	
	# load the annotation
	dannot, dgid2tid = load_info(args.info)

	# update the kpseudo with target and gene-level downweigthing
	message("Applying gene ambiguity down-weighting to equivalence classes")
	kp.set_targets(tid2faiIndex, lfai)
	kp.calc_ec_weights(dannot)
	kp.calc_sample_target_counts()

	# load the sample abundance files
	sample_abundance = [None for i in range(kp.num_samples)]
	i = 0
	for sid in kp.samples:

		if os.path.isdir(sid):
			message("Found sample {}".format(sid))
			sample_abundance[i] = KAbundance("{}/abundance.tsv".format(sid))

		i += 1


	#
	# now using the dgid2tid table we can gather the relative isoform levels from the 
	# abundance files and the total count from the pseudo quantification. the relative
	# levels are used to modify the psuedo based count for each isoform.
	#
	
	message("Using relative isoform count abundances to augment the pseudo counts")
	for gid in dgid2tid.keys():
		num_tid = len(dgid2tid[gid])
		if num_tid > 1:
			# visit each sample, gather the relative isoform levels and the total counts 
			# from the pseudo output
			for j in range(kp.num_samples):
				iso_vals = []
				total_count = 0
				for tid in dgid2tid[gid]:
					# get index from fai
					fai_index = tid2faiIndex[tid]
					# get count
					total_count += kp.sample_target_counts[j][fai_index]
					# get kallisto estimate
					iso_vals.append(sample_abundance[j].get_count(tid))

				sum_count = sum(iso_vals)
				if sum_count > 0:
					for i in range(num_tid):
						iso_vals[i] = iso_vals[i]/sum_count * total_count

				# now plug the counts back into pseudo target table
				for i in range(num_tid):
					tid = dgid2tid[gid][i]
					faidx = tid2faiIndex[tid]
					kp.sample_target_counts[j][faidx] = iso_vals[i]


	# 
	# build a final output file per sample
	#
	message("writing per-sample *.kpquant files")
	# tid, gid, gname, locus, length, eff_length, est_hits
	for i in range(kp.num_samples):

		with open("{}.kpquant".format(kp.samples[i]), "w") as fout:
			fout.write("transcript_id\tgene_id\tgene_name\tlocus\tlength\teff_length\test_count\n")

			# use the fai as a guide
			for j in range(len(lfai)):
				tid = lfai[j]

				gene_id = dannot[tid]['gene_id']
				gene_name = dannot[tid]['gene_name']
				locus = dannot[tid]['locus']
				length = sample_abundance[i].get_length(tid)
				eff_length = sample_abundance[i].get_eff_length(tid)
				count = kp.sample_target_counts[i][tid2faiIndex[tid]]
				est_count = 0
				if eff_length > 0:
					est_count = count * length/eff_length


				lout = [tid, gene_id, gene_name, locus, length, eff_length, est_count]
				fout.write("\t".join(map(str, lout)))
				fout.write("\n")



	return 0


def load_fai(f):

	lfai = []
	tid2index = {}
	index = 0

	with open(f, "r") as fin:

		for szl in fin:
			aln = szl.strip().split("\t")
			tid = aln[0]
			lfai.append(tid)
			tid2index[tid] = index
			index += 1

	message("Loaded FAI. {} targets".format(index))

	return lfai, tid2index

def load_info(f):

	dtid = {}
	gid2tid = {}
	count = 0

	with open(f, "r") as fin:

		for szl in fin:
			aln = szl.strip().split("\t")
			if aln[0] == "chrom":
				# header!
				continue

			count += 1
			tid = aln[3]
			gid = aln[2]

			dtid[tid] = { "chrom": aln[0], "gene_id": aln[2], "gene_name": aln[4], "locus": aln[5] }

			if gid not in gid2tid:
				gid2tid[gid] = []

			gid2tid[gid].append(tid)

	message("Loaded annotation. {} targets".format(count))

	return dtid, gid2tid


def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

#
# klass for kallisto pseudo stuff
class KPseudo(object):

	def __init__(self, path):
		self.path = path
		self.ready = False

		self.samples = []
		self.num_samples = 0
		self.ec = []
		self.ec_weights = []
		self.num_ec = 0
		self.num_targets = 0
		self.sample_ec_counts = []
		self.sample_target_counts = []
		self.sample_totals = []
		self.target_map = None
		self.targets = None
		
	def check_path(self):
		# check the path for the necessary files
		
		if self.ready:
			return True

		dlist = set(os.listdir(self.path))
		ready = True
		expected = ['matrix.cells', 'matrix.ec', 'matrix.tsv']
		for f in expected:
			if f not in dlist:
				ready = False

		self.ready = ready
		return ready

	def load_data(self):

		# load the samples up
		with open("{}/matrix.cells".format(self.path), "r") as fin:
			for szl in fin:
				szl = szl.strip()
				self.samples.append(szl)

		message("found {} samples".format(len(self.samples)))
		self.num_samples = len(self.samples)
		self.sample_totals = [0 for i in range(self.num_samples)]

		with open("{}/matrix.ec".format(self.path), "r") as fin:
			for szl in fin:
				aln = szl.strip().split("\t")
				self.ec.append(map(int, aln[1].split(",")))
				tmp = max(self.ec[-1])+1
				if tmp > self.num_targets:
					self.num_targets = tmp

		self.num_ec = len(self.ec)
		# ec weights used when creating target level counts
		self.ec_weights = [0 for i in range(self.num_ec)]

		message("found {} equivalence classes and {} targets".format(len(self.ec), self.num_targets))

		self.sample_ec_counts = [[] for i in range(self.num_samples)]
		self.sample_target_counts = [[] for i in range(self.num_samples)]
		
		for i in range(self.num_samples):
			# pre-populate with zeros
			self.sample_ec_counts[i] = [0 for j in range(self.num_ec)]
			self.sample_target_counts[i] = [0.0 for j in range(self.num_targets)]

		# now we can read in the 'tsv' file
		with open("{}/matrix.tsv".format(self.path), "r") as fin:
			for szl in fin:
				aln = map(int, szl.strip().split("\t"))
				ecid = aln[0]
				# update sample total count
				self.sample_totals[aln[1]] += aln[2]
				self.sample_ec_counts[aln[1]][ecid] = aln[2]
				wi = 1.0/len(self.ec[ecid])
				self.ec_weights[ecid] = wi
				for tid in self.ec[ecid]:
					self.sample_target_counts[aln[1]][tid] += float(aln[2])*wi

	def set_targets(self, d, l):
		self.target_map = d
		self.targets = l
		return 0
	
	def calc_ec_weights(self, annot):
		if self.target_map is None:
			message("You have to set the target dict to translate transcript ids to target indices first")
			return 1
		
		# loop through the EC. at each EC gather not only the count of transcripts but also 
		# the number of genes.
		
		# initalize the ec weight vector
		self.ec_weights = [0 for i in range(self.num_ec)]
		for i in range(self.num_ec):
			num_targets = len(self.ec[i])
			# get gene names
			gnames = set()
			for tidx in self.ec[i]:
				tid = self.targets[tidx]
				gnames.add(annot[tid]['gene_name'])
			
			gnames = list(gnames)
			num_genes = len(gnames)
			# calculate weight
			self.ec_weights[i] = 1.0/num_targets * 1.0/(num_genes**2)
		
		return 0
			
	def calc_sample_target_counts(self):
		
		for i in range(self.num_samples):
			# pre-populate with zeros
			self.sample_target_counts[i] = [0.0 for j in range(self.num_targets)]
		
		for i in range(self.num_ec):
			tidx = self.ec[i]
			wi = self.ec_weights[i]
			for j in range(self.num_samples):
				for k in tidx:
					self.sample_target_counts[j][k] += float(self.sample_ec_counts[j][i])*wi
		
		return 0

	def get_total_sample_counts(self):
		tmp = [0 for i in range(self.num_samples)]
		for i in range(self.num_samples):
			for j in range(self.num_ec):
				tmp[i] += self.sample_ec_counts[i][j]

		return tmp

class KAbundance(object):

	def __init__(self, f):
		
		self.fname = f
		self.d = {}

		self.lengths = []
		self.eff_lengths = []
		self.counts = []
		self.est_counts = []

		self.load()


	def load(self):

		index = 0
		message("Loading {}".format(self.fname))
		with open(self.fname, "r") as fin:
			# skip header: tid, length, eff_length, est_counts, tpm
			szl = fin.readline()
			for szl in fin:
				aln = szl.strip().split("\t")

				tid = aln[0]
				length = float(aln[1])
				eff_length = float(aln[2])
				est_count = float(aln[3])
				count = 0
				if est_count > 0:
					# un-effective the count
					count = est_count * eff_length/length

				self.d[tid] = index
				self.lengths.append(length)
				self.eff_lengths.append(eff_length)
				self.counts.append(count)
				self.est_counts.append(est_count)
				index += 1

		return 0

	def get_count(self, tid):
		return self.counts[self.d[tid]]

	def get_length(self, tid):
		return self.lengths[self.d[tid]]

	def get_eff_length(self, tid):
		return self.eff_lengths[self.d[tid]]



def print_exception():
	exc_type, exc_value, exc_traceback = sys.exc_info()
	print "*** print_tb:"
	traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
	print "*** print_exception:"
	traceback.print_exception(exc_type, exc_value, exc_traceback,
	                          limit=2, file=sys.stdout)
	print "*** print_exc:"
	traceback.print_exc()
	print "*** format_exc, first and last line:"
	formatted_lines = traceback.format_exc().splitlines()
	print formatted_lines[0]
	print formatted_lines[-1]
	print "*** format_exception:"
	print repr(traceback.format_exception(exc_type, exc_value,
	                                      exc_traceback))
	print "*** extract_tb:"
	print repr(traceback.extract_tb(exc_traceback))
	print "*** format_tb:"
	print repr(traceback.format_tb(exc_traceback))
	print "*** tb_lineno:", exc_traceback.tb_lineno	

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('fai', type=str, help="FAI index of the FASTA used to build the kallisto index used for quantification.")
parser.add_argument('info', type=str, help="'info' file with per-transcript annotation in the usual format.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

