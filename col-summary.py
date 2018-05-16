#!/usr/bin/python
#==============================================================================
# col-summary.py
#
# Shawn Driscoll
# 20170622
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Produce column summaries of a tab delim data file
#==============================================================================

import sys
import argparse
import math
import re
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
	
	ncol = 0
	
	ldata = []
	lheader = []
	
	if args.file == "-":
		fin = sys.stdin
	else:
		fin = open(args.file, "r")
	
	# parse the file into lists by columns
	if args.s > 0:
		i = 0
		while i < args.s:
			szl = fin.readline()
			if i==0:
				lheader = szl.strip().split("\t")
			i += 1
	
	for szl in fin:
		aln = szl.strip().split("\t")
		if len(ldata)==0:
			ncol = len(aln)
			ldata = [[] for i in range(ncol)]
		
		for j in range(ncol):
			ldata[j].append(aln[j])
			
	# close file
	fin.close()			

	if len(lheader) == 0:
		for i in range(ncol):
			lheader.append("col{}".format(i+1))

	ss = Stats()
	print "\t".join(ss.tolist(header=True))

	for i in range(ncol):
		ss = Stats()
		ss.collect(ldata[i])
		ss.calculate()
		tmp = [lheader[i]] + ss.tolist()
		print "\t".join(map(str, tmp))	

	return 0

def summarize(v):
	n = len(v)
	
class Stats(object):
	def __init__(self):
		self.n = 0
		
		self.levels = None
		self.values = None
		
		self.n_values = 0
		self.n_levels = 0
		
		self.mean = 0
		self.stdev = 0
		self.cv = 0
		
		self.median = 0
		self.q25 = 0
		self.q75 = 0
		self.iqr = 0
		
		self.sum = 0
		self.max = 0
		self.min = 0
		
		self.unique_levels = 0
		self.max_level_count = 0
		self.max_level = "NA"
		self.min_level_count = 0
		self.min_level = "NA"		
	
	# calculate stats on vector, v
	def collect(self, v):
		self.n = len(v)
		
		self.values = []
		self.levels = defaultdict(int)
		
		for i in range(self.n):
			x = v[i]
			try:
				x = float(x)
				# we have a number!
				self.values.append(x)
				self.n_values += 1
			except:
				# not a number
				self.levels[x] += 1
				self.n_levels += 1
		
		return 0
		
	def calculate(self):
		
		if self.n_values > 0:
			
			v = np.array(self.values)
			self.sum = np.sum(v)
			self.mean = np.mean(v)
			self.stdev = np.std(v)
			if self.mean != 0:
				self.cv = self.stdev/self.mean
			
			quants = np.percentile(v, [0.25, 0.5, 0.75])
			self.median = quants[1]
			self.q25 = quants[0]
			self.q75 = quants[2]
			self.iqr = quants[2]-quants[1]
			
			self.min = np.min(v)
			self.max = np.max(v)
		
		if self.n_levels > 0:
			kk = self.levels.keys()
			self.unique_levels = len(kk)
			
			# initalize
			self.min_level_count = self.levels[kk[0]]
			self.min_level = kk[0]
			
			for k in kk:
				if self.levels[k] > self.max_level_count:
					self.max_level = k
					self.max_level_count = self.levels[k]
				
		
		return 0
			
	def tolist(self, header=False):
		
		if header:
			lout = ["N", "N_values", "N_factors", "sum", "mean", "stdev", "cv", "median", "q25", "q75", "iqr", "min", "max", "U_unique_factors", "max_factor_count", "min_factor_count", "max_factor", "min_factor"]
		
		else:
				
			lout = [self.n, self.n_values, self.n_levels]
			if self.n_values > 0:
				lout += [self.sum, self.mean, self.stdev, self.cv, self.median, self.q25, self.q75, self.iqr, self.min, self.max]
			else:
				lout += ["NA" for i in range(10)]
			
			if self.n_levels > 0:
				lout += [self.unique_levels, self.max_level_count, self.min_level_count, self.max_level, self.min_level]
			else:
				lout += ["NA" for i in range(5)]
		
		return lout
				
		

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


parser = argparse.ArgumentParser(description="Produce column summaries of a tab delim data file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file', type=str, help="Input file. Must be tab delim")
parser.add_argument('-s', type=int, default=0, help="Skip this number of lines from top of file.")
parser.add_argument('-z', action="store_const", const=True, default=False, 
	help="Ignore zeros in numeric summaries")
	
args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

