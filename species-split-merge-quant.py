#!/usr/bin/python
#==============================================================================
# species-split-merge-quant.py
#
# Shawn Driscoll
# 20170713
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Merges 'bam-gcounts-cs' type quantifications of the two species from a 
# split experiment. The expressions for the second species are going to be 
# made relative to the first. This way if we're dealing with a transplanted 
# species we can see the transplanted expressions relative to the host across
# samples making it possible to determin the relative amounts of transplanted
# RNA per sample. 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
from time import localtime
from collections import defaultdict
import copy

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

OUT_TID = 0
OUT_GID = 1
OUT_GNAME = 2
OUT_CHROM = 3
OUT_STRAND = 4
OUT_LENGTH = 5
OUT_EFF_LENGTH = 6
OUT_HITS = 7
OUT_EFF_HITS = 8
OUT_FPKM = 9
OUT_TPM = 10

#==============================================================================
# main
#==============================================================================


def main(args):
	
	#
	# variables
	#
	
	quantA = None
	quantA0 = None
	
	#
	# load the two quants
	#
	
	if not isfile(args.quantA):
		message("Input file does not exist {}".format(args.quantA))
		return 1
	if not isfile(args.quantB):
		message("Input file does not exist {}".format(args.quantB))
		return 1
	
	quantA = QuantFile()
	with open(args.quantA, "r") as fin:
		# skip header
		sz = fin.readline()
		
		quantA.add_header(sz)
		
		for sz in fin:
			qline = QuantRow(sz)
			quantA.add_line(qline)
		
	# continue by adding lines from the second file
	with open(args.quantB, "r") as fin:
		# header
		sz = fin.readline()
		
		for sz in fin:
			qline = QuantRow(sz)
			quantA.add_line(qline, exclude_from_totals=True)
	
	# update FPKM and TPM 
	quantA.update_fpkm()
	quantA.update_tpm()
	
	# write results
	print quantA.header
	for qrow in quantA.lines:
		print "\t".join(qrow.to_list())
	
	return 0


#==============================================================================
# general functions
#==============================================================================


def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

class QuantRow(object):
	def __init__(self, sz):
		aln = sz.strip().split("\t")
		
		self.transcript_id = aln[0]
		self.gene_id = aln[1]
		self.gene_name = aln[2]
		self.chrom = aln[3]
		self.strand = aln[4]
		self.length = float(aln[5])
		self.eff_length = float(aln[6])
		self.hits = float(aln[7])
		self.eff_hits = float(aln[8])
		self.fpkm = float(aln[9])
		self.tpm = float(aln[10])
		
	def __str__(self):
		return "\t".join(self.to_list())

	def to_list(self):
		lout = [self.transcript_id, self.gene_id, self.gene_name, self.chrom, self.strand, 
			str(self.length), "{:0.4f}".format(self.eff_length), "{:0.4f}".format(self.hits), 
			"{:0.4f}".format(self.eff_hits), "{:0.4f}".format(self.fpkm), "{:0.4f}".format(self.tpm)]
		return lout

class QuantFile(object):
	
	def __init__(self):
		self.header = ""
		self.lines = []
		
		self.num_lines = 0
		self.total_hits = 0
		self.total_eff_hits = 0
		self.total_fpkm = 0
		
		self.exclude = []

	def add_header(self, sz):
		self.header = sz.strip()
		return 0
	
	# append a quant object line to the file
	def add_line(self, q, exclude_from_totals=False):
		self.lines.append(q)
		self.num_lines += 1
		
		self.exclude.append(exclude_from_totals)
		
		if not exclude_from_totals:
			self.total_hits += q.hits
			self.total_eff_hits += q.eff_hits
			self.total_fpkm += q.fpkm
		
		return 0
	
	def update_total_hits(self):
		self.total_hits = 0
		for i in range(self.num_lines):
			if not self.exclude[i]:
				self.total_hits += self.lines[i].hits
		return 0
	
	def update_eff_hits(self):
		self.total_eff_hits = 0
		for i in range(self.num_lines):
			if self.lines[i].eff_length > 0:
				self.lines[i].eff_hits = self.lines[i].hits*float(self.lines[i].length)/self.lines[i].eff_length
				# update total?
				if not self.exclude[i]:
					self.total_eff_hits += self.lines[i].eff_hits
			else:
				self.lines[i].eff_hits = 0
		return 0
		
	def update_fpkm(self):
		self.total_fpkm = 0
		for i in range(self.num_lines):
			f = self.lines[i].eff_hits*1e9/(float(self.lines[i].length)*self.total_eff_hits)
			# update total
			if not self.exclude[i]:
				self.total_fpkm += f
			self.lines[i].fpkm = f
		return 0
	
	def update_tpm(self):
		for i in range(self.num_lines):
			self.lines[i].tpm = self.lines[i].fpkm*1e6/self.total_fpkm
		return 0
		

def header_to_dict(sz):
	aln = sz.strip().split("\t")
	d = defaultdict(int)
	i = 0
	for k in aln:
		d[k] = i
		i += 1
	
	return d

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Combine split-species quantifications making TPM levels that are relative to the FIRST of the two files.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
parser.add_argument('quantA', type=str, help="First or host quantification")
parser.add_argument('quantB', type=str, help="Second or guest quantification")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


