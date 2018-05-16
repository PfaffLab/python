#!/usr/bin/python
#==============================================================================
# gtf-assembly-pipeline.py
#
# Shawn Driscoll
# 20160427
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# given a TSV file with two columns: assembly name (gtf) \t  condition
# this script runs gffcompare and stringtie --merge to generate a final
# version.
#==============================================================================

import sys, argparse 
import math 
import re
import copy
from os.path import isfile, expanduser
import subprocess as sp
import os

# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
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
	conditions = {}
	filter_tracking = False
	condition_gtfs = []

	# check options
	if args.f != 0:
		filter_tracking = True

	# --
	# load the samples from the tsv file
	fin = open(args.gtf_list, "r")
	szl = fin.readline()
	for szl in fin:
		aln = szl.strip().split("\t")
		
		# add sample into the conditions dict
		if aln[1] not in conditions:
			conditions[aln[1]] = []

		conditions[aln[1]].append(aln[0])

	fin.close()

	# --
	# run gffcompare per condition
	# --

	for cond in conditions.keys():
		run_gffcompare(conditions[cond], args.ref, cond)

		# filter tracking
		if filter_tracking:
			sys.stderr.write("filtering {} tracking file...\n".format(cond))
			# need to filter the tracking file and then the GTF file as well so a dict
			# of transcript ids to keep will be made
			diso = {}
			fin = open("{}.tracking".format(cond), "r")
			trow = TrackingRow()
			for szl in fin:
				trow.__init__()
				trow.parse(szl.strip())

				if args.f > 0 and args.f < 1:
					# use ratio of missing
					mrat = 1 - trow.present_samples*1.0/trow.num_samples
					if mrat < args.f:
						# keep it
						diso[trow.tid] = 0
				elif args.f >= 1:
					# max number of samples missing
					mmis = trow.num_samples - trow.present_samples
					if mmis <= args.f:
						# keep it
						diso[trow.tid] = 0

			fin.close()

			sys.stderr.write("applying filtering to {} combined GTF...\n".format(cond))

			# open GTF, print back out only the rows that pass the filter
			fin = open("{}.combined.gtf".format(cond), "r")
			fout = open("{}.combined.filtered.gtf".format(cond), "w")
			grow = GtfRow()
			for szl in fin:
				grow.__init__()
				grow.parse(szl.strip())
				if grow.tid in diso:
					# print this back out
					fout.write(szl)

			fin.close()
			fout.close()

			condition_gtfs.append("{}.combined.filtered.gtf".format(cond))
		else:
			condition_gtfs.append("{}.combined.gtf".format(cond))


	# --
	# run stringtie merge
	# --

	merge_stub = "_".join(conditions.keys())
	run_merge(condition_gtfs, args.ref, merge_stub)

	# --
	# run the final gffcompare
	# --
	
	run_gffcompare(["{}.gtf".format(merge_stub)], args.ref, args.stub)

	# --
	# clean up!
	# --


	for cond in conditions.keys():
		os.unlink("{}.tracking".format(cond))
		os.unlink("{}.loci".format(cond))
		os.unlink("{}.stats".format(cond))

	os.unlink("{}.gtf".format(merge_stub))


	return 0


#==============================================================================
# classes
#==============================================================================

# class for *.tracking file rows
class TrackingRow(object):
	def __init__(self):
		self.tid = ""
		self.xloc = ""
		self.ref_gene = ""
		self.ref_tid = ""
		self.class_code = ""
		self.samples = []
		self.num_samples = 0
		self.present_samples = 0
		return None

	# 
	# this function parses the tracking file row into this class
	def parse(self, sz):
		
		aln = sz.strip().split("\t")
		n = len(aln)

		if aln[2]=="-":
			ref_split = ["-", "-"]
		else:
			ref_split = aln[2].split("|")

		self.tid = aln[0]
		self.xloc = aln[1]
		self.ref_gene = ref_split[0]
		self.ref_tid = ref_split[1]
		self.class_code = aln[3]

		# copy the samples into a list
		self.samples = list(aln[4:n])
		# get number of samples
		self.num_samples = len(self.samples)
		# get count of samples that contain the feature
		for i in range(self.num_samples):
			if self.samples[i] != "-":
				self.present_samples += 1

		return 0

	def tolist(self):
		lout = [
			self.tid, 
			self.xloc, 
			"|".join([self.ref_gene, self.ref_tid]), 
			self.class_code
		]
		lout += self.samples
		return lout


#
# class to hold GTF row objects. 
class GtfRow(object):
	def __init__(self, rname="", db="", type="", start=0, end=0, strand="."):
		self.rname = rname
		self.db = db
		self.type = type
		self.start = start
		self.end = end
		self.strand = strand
		
		self.tid = ""
		self.gid = ""
		
		self.attrs = {}

		return None
	
	# parse info in GTF row, sz, into this object
	def parse(self, sz):
		aln = sz.strip().split("\t")
		self.rname = aln[0]
		self.db = aln[1]
		self.type = aln[2]
		self.start = int(aln[3])
		self.end = int(aln[4])
		self.strand = aln[6]
		
		# parse attributes from field 8
		fsplit = aln[8].split("\"")
		n = len(fsplit)-1
		i = 0
		while i < n:
			key = re.sub(';','',fsplit[i])
			self.attrs[key.strip()] = fsplit[i+1].strip()
			i += 2

		self.tid = self.attrs['transcript_id']
		self.gid = self.attrs['gene_id']
		
		return 0

	#--
	# tolist
	# return a list version of this with elements in place of the columns 
	# of a GTF
	def tolist(self):
		ltmp = [self.rname, self.db, self.type, 
			int(self.start), int(self.end), ".", self.strand, "."]
		
		szattr = "transcript_id \"{}\"; gene_id \"{}\";".format(self.tid, self.gid)
		akey = self.attrs.keys()
		if len(akey) > 0:
			# append additional attributes to the szattr string
			for aid in akey:
				szattr += " {} \"{}\";".format(aid, self.attrs[aid])
		
		ltmp.append(szattr)
		
		return(ltmp)

#==============================================================================
# functions
#==============================================================================

def run_gffcompare(gtf_list, ref, stub):
	cmd = "gffcompare -r {} -TC -o {} {}".format(ref, stub, " ".join(gtf_list))
	return runcmd(cmd)

def run_merge(gtf_list, ref, stub):
	cmd = "stringtie --merge -G {} -T 0 -F 0 {}".format(ref, " ".join(gtf_list))
	
	sys.stderr.write("CMD: {} > {}.gtf\n".format(cmd, stub))

	fout = open("{}.gtf".format(stub), "w")

	p1 = sp.Popen(cmd.split(), stdout=fout)
	p1.wait()
	return(0)


# --
# runcmd
# run a system level command in subprocess. optionally you can return the process.
# if the process isn't returned then the function waits for the process to finish
def runcmd(cmd, returnProcess=False):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="GTF assembly pipeline")
parser.add_argument('stub', type=str, help="Name stub for the build")
parser.add_argument('gtf_list', type=str, 
	help="Tab-delim file with two columns for the gtf name and condition. expected to have a header line.")
#parser.add_argument('ref', type=str, nargs="+", help="Junctions file(s)")
parser.add_argument("ref", type=str, 
				help="Known annotation (gtf)")
parser.add_argument("-f", default=0, type=float, 
	help="Tracking file filtering. Max # samples without transcript (int) or max ratio of missing (float)")


args = parser.parse_args()

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

