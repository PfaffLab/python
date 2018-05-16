#!/usr/bin/python
#==============================================================================
# species-split-combine-stats.py
#
# Shawn Driscoll
# 20170615
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Combine stats from multiple species-split-ab.py
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser

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
	n = len(args.runlogs)
	fout_stats = "species_split_combined_stats.tsv"
	
	logs = MultiLog()
	
	# load the logs
	for f in args.runlogs:
		log = SplitLog()
		log.parse(f)	
		logs.add_log(log)
	
	with open(fout_stats, "w") as fout:
		for l in logs.get_stats_table():
			fout.write("\t".join(l) + "\n")	
	
	return 0

class SplitLog(object):
	
	def __init__(self):
		self.file = None
		self.sample = None
		self.stat_order = []
		self.stats = {}
		
	
	def parse(self, fname):
		self.file = fname
		
		with open(fname, "r") as fin:
			# skip header
			szl = fin.readline()
			for szl in fin:
				aln = szl.strip().split("\t")
				self.sample = aln[0]
				self.stat_order.append(aln[1])
				self.stats[aln[1]] = int(aln[2])
		
		return 
			
			
class MultiLog(object):
	def __init__(self):
		self.logs = []
		self.stat_order = []
		self.num_logs = 0
	
	def add_log(self, obj):
		self.logs.append(obj)
		self.num_logs += 1
		
		if len(self.stat_order)==0:
			self.stat_order = list(obj.stat_order)
	
	def get_stats_table(self):
		
		if self.num_logs==0:
			return None
		
		outlines = []
		# append header
		lout = ["sample"] + self.stat_order
		outlines.append(lout)
		
		# loop through logs
		for log in self.logs:
			lout = [log.sample]
			for stat in self.stat_order:
				lout.append(log.stats[stat])
			
			outlines.append(map(str, lout))
		
		return outlines
				

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Combine stats from multiple 'species-split-ab.py' outputs.")
parser.add_argument('runlogs', type=str, nargs="+", 
	help="Specify 2 or more 'tsv' files from 'species-split-ab.py' output")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

