#!/usr/bin/python
#==============================================================================
# scrna-bd-combine-results.py
#
# Shawn Driscoll
# 20170615
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Combine quantification results from multiple scrna-gcounts-bd runs
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime

# from igraph import *
# from subprocess import Popen
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
	n = len(args.schits)
	allQ = QuantCatalog()
	has_bulk = False
	
	stub = "scrna_"

	if args.bulk is not None:
		message("Loading bulk expression data")
		blk = BulkSample()
		blk.parse(args.bulk)
		message("Creating quantiles")
		blk.calculate_quantiles()
		has_bulk = True
			
	message("Loading files")
	for f in args.schits:
		progress_message("{}".format(f))
		quant = scRnaHits()
		quant.parse(f)
		
		if has_bulk:
			# evaluate sensitivity
			quant.calc_sensitivity(blk)
				
		allQ.add_quant(quant)
	
	sys.stderr.write("\n")
	
	
	# all quants are loaded. now make some stats and generate the quantificaiton 
	# summary tables
	message("building tables")
	allQ.build_tables()
	
	#
	# write tables
	#
	
	message("Writing raw reads table")
	with open("{}raw_reads.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.read_table)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.read_table[i])))
			fout.write("\n")

	message("Writing raw MI table")
	with open("{}raw_MI.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.read_table)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.raw_mi_table[i])))
			fout.write("\n")

	message("Writing MI table")
	with open("{}corrected_MI.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.read_table)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.mi_table[i])))
			fout.write("\n")
	
	message("Writing summary report")
	
	rpt = allQ.build_report()
	with open("{}report.tsv".format(stub), "w") as fout:
		for l in rpt:
			fout.write("\t".join(map(str, l)) + "\n")

	return 0

class QuantCatalog(object):
	def __init__(self):
		self.quants = []
		self.num_quants = 0
	
		self.annot = None
		self.annot_header = ["gene_id", "gene_name", "chrom", "strand"]
		self.read_table = None
		self.raw_mi_table = None
		self.mi_table = None
		self.data_header = None
		
		
	def add_quant(self, obj):
		self.num_quants += 1
		self.quants.append(obj)
	
	def build_tables(self):
		
		# data header is all of the sample names
		self.data_header = [obj.sample for obj in self.quants]
				
		# build a list of all gene id from all samples
		all_gid = set()
		for obj in self.quants:
			all_gid.update(obj.gid2idx.keys())
		all_gid = sorted(list(all_gid))
		
		num_features = len(all_gid)
		# create tables of 'num_features' lists of length 'self.num_quants'
		self.annot = [[] for j in range(num_features)]
		self.read_table = [[0 for i in range(self.num_quants)] for j in range(num_features)]
		self.raw_mi_table = [[0 for i in range(self.num_quants)] for j in range(num_features)]
		self.mi_table = [[0 for i in range(self.num_quants)] for j in range(num_features)]

		# loop through gene id and populate all of the tables
		
		for j in range(len(all_gid)):
			gid = all_gid[j]
			for i in range(self.num_quants):
				obj = self.quants[i]
				if gid in obj.gid2idx:
					idx = obj.gid2idx[gid]
					self.annot[j] = obj.annot[idx]
					self.read_table[j][i] = obj.reads[idx]
					self.raw_mi_table[j][i] = obj.raw_mi[idx]
					self.mi_table[j][i] = obj.mi[idx]
					
		return
		
	def build_report(self):
		
		# figure out if the quants have bulk data associations
		has_bulk = True
		for obj in self.quants:
			if not obj.has_sense:
				has_bulk = False
		
		# build initial table
		rtable = [[] for j in range(self.num_quants+1)]
		
		rtable[0] = ["Sample", "Reads", "Raw_MI", "MI", "Genes_Over_0_MI", "Genes_Over_1_MI"]
		if self.quants[0].has_sense:
			has_sense = True
		
		if has_bulk:
			rtable[0] += ["q0", "q10", "q20", "q30", "q40", "q50", "q60", "q70", "q80", "q90"]
		
		for i in range(self.num_quants):
			obj = self.quants[i]
			j = i+1
			rtable[j] = [obj.sample, obj.total_reads, obj.total_raw_mi, obj.total_mi, obj.mi_over1, obj.mi_over2]
			if has_bulk:
				for k in range(len(obj.bulk_bin_exp)):
					# rtable[j].append(obj.bulk_bin_exp[k]*1.0/obj.mi_over1)
					sz = "{}:{}:{}".format(obj.bulk_bin_exp[k], obj.bulk_bin_pres[k], obj.bulk_bin_size[k])
					rtable[j].append(sz)
					
		
		return rtable

class scRnaHits(object):
	
	def __init__(self):
		self.file = None
		self.sample = None
		
		self.gid2idx = {}
		self.gn2idx = {}
		self.annot = []
		self.reads = []
		self.raw_mi = []
		self.mi = []
		
		self.total_reads = 0
		self.total_raw_mi = 0
		self.total_mi = 0
		
		self.mi_over1 = 0
		self.mi_over2 = 0
		
		self.has_sense = False
		self.bulk_bin_size = None
		self.bulk_bin_exp = None
		self.bulk_bin_pres = None
	
	def parse(self, fname):
		self.file = fname
		
		tmp = fname.split("/")
		tmp2 = tmp[-1].split(".")
		self.sample = tmp2[0]
		
		# parse the file
		idx = 0
		with open(fname, "r") as fin:
			# skip header
			szl = fin.readline()
			# parse on
			for szl in fin:
				aln = szl.strip().split("\t")
				# note gene id index
				self.gid2idx[aln[0]] = idx
				if aln[1] not in self.gn2idx:
					self.gn2idx[aln[1]] = []
				
				self.gn2idx[aln[1]].append(idx)
				
				idx += 1
				
				self.annot.append(list(aln[0:4]))
				self.reads.append(int(aln[4]))
				self.raw_mi.append(int(aln[5]))
				self.mi.append(int(aln[6]))
				
				if int(aln[6]) > 0:
					self.mi_over1 += 1
				if int(aln[6]) > 1:
					self.mi_over2 += 1
		
		self.total_reads = sum(self.reads)
		self.total_raw_mi = sum(self.raw_mi)
		self.total_mi = sum(self.mi)
		
		return
	
	def calc_sensitivity(self, obj):
		# loop through the quantile bins in obj and count number of genes in each
		# bin that appear to be expressed in this sample
		
		num_bins = len(obj.quantile_bins)
		
		bin_exp = [0 for i in range(num_bins)]
		bin_present = [0 for i in range(num_bins)]
		# copy the bin sizes
		bin_size = [obj.bin_size[i] for i in range(num_bins)]
		
		for i in range(len(obj.quantile_bins)):
			for gname in obj.quantile_bins[i]:
				if gname in self.gn2idx:
					idx = self.gn2idx[gname]
					bin_present[i] += 1
					for j in idx:
						if self.mi[j] > 0:
							bin_exp[i] += 1
							# done with this gene
							break
		
		self.bulk_bin_pres = bin_present
		self.bulk_bin_exp = bin_exp
		self.bulk_bin_size = bin_size
		self.has_sense = True
		
		return

class BulkSample(object):
	
	def __init__(self):
		
		self.id2idx = {}
		
		self.id = []
		self.mean_tpm = []
		self.quantiles = []
		self.quantile_bins = []
		self.bin_size = []
	
	def parse(self, fname):
		
		idx = 0
		with open(fname, "r") as fin:
			# skip header
			szl = fin.readline()
			for szl in fin:
				aln = szl.strip().split("\t")
				
				if aln[0] not in self.id2idx:
					self.id2idx[aln[0]] = idx
					self.id.append(aln[0])
					
					if len(aln) > 2:
						tpm = map(float, aln[1:len(aln)])
						self.mean_tpm.append(sum(tpm)/len(tpm))
					else:
						self.mean_tpm.append(float(aln[1]))
				else:
					idx = self.id2idx[aln[0]]
					if len(aln) > 2:
						tpm = map(float, aln[1:len(aln)])
						self.mean_tpm[idx] += sum(tpm)/len(tpm)
					else:
						self.mean_tpm[idx] += float(aln[1])
	
	def calculate_quantiles(self):
		# build a table of non-zero tpm with their indices
		ltpm = []
		idx = 0
		
		for i in range(len(self.mean_tpm)):
			if self.mean_tpm[i] > 0.001:
				ltpm.append([i, self.mean_tpm[i]])
		
		# get total rows
		n = len(ltpm)
		
		# sort this thing by the tpm values
		ltpm.sort(key=lambda x:x[1])
		
		# start the quantile vector
		self.quantiles = [0 for i in range(len(self.mean_tpm))]
		
		# loop through and assign the quantile index
		for i in range(n):
			self.quantiles[ltpm[i][0]] = (i+1)*1.0/(n+1)
		
		# now with the quantiles we can bin them into 10's
		
		self.quantile_bins = [[] for i in range(10)]
		self.bin_size = [0 for i in range(10)]
		
		for i in range(len(self.quantiles)):
			# only assign non-zero
			if self.mean_tpm[i] > 0.001:
				qadj = int(math.floor(self.quantiles[i]*10))
				self.quantile_bins[qadj].append(self.id[i])
				self.bin_size[qadj] += 1
		
		return

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Combine multiple scRNA quantifications from 'scrna-gcounts-bd.py'")
parser.add_argument('schits', type=str, nargs="+", 
	help="Specify 2 or more 'schits' files from 'scrna-gcounts-bd.py' output")

parser.add_argument('-b', '--bulk', type=str, default=None, 
	help="Bulk RNA-Seq gene level TPM quantification to evaluate sensitivity of the scRNA data")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

