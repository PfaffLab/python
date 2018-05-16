#!/usr/bin/python
#==============================================================================
# scrna-bd-combine-stats.py
#
# Shawn Driscoll
# 20170615
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Combine read alignment/assignment stats from multiple 'scrna-gcounts-bd.py' 
# runs.  
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
	n = len(args.runlogs)
	fout_stats = "scrna_combined_stats.tsv"
	fout_mapq = "scrna_combined_mapq.tsv"
	
	logs = MultiRunLog()
	
	# load the logs
	for f in args.runlogs:
		log = RunLog()
		log.parse(f)	
		logs.add_log(log)
	
	with open(fout_stats, "w") as fout:
		for l in logs.get_stats_table():
			fout.write("\t".join(l) + "\n")	
	
	with open(fout_mapq, "w") as fout:
		for l in logs.get_mapq_table():
			fout.write("\t".join(l) + "\n")

	return 0

class  RunLog(object):
	
	def __init__(self):
		self.file = None
		self.bam = None
		self.runtime = None
		self.stat_order = []
		self.stats = {}
		self.mapq = {}
	
	def parse(self, fname):
		self.file = fname
		
		with open(fname, "r") as fin:
			# read in the first line - it's nothing
			szl = fin.readline()
			# read in run time
			self.runtime = re.sub("^# ", "", fin.readline().strip())
			# read in bam name
			self.bam = re.sub("^# ", "", fin.readline().strip())
			# read 3 lines
			skip = 3
			while skip > 0:
				szl = fin.readline()
				skip -= 1
			
			# loop through stats
			while True:
				szl = fin.readline().strip()
				if len(szl)==0 or szl=="":
					break
					
				aln = szl.split("\t")
				self.stat_order.append(aln[0])
				self.stats[aln[0]] = int(aln[1])
				
			# skip to the MAPQ table
			skip = 2
			while skip > 0:
				szl = fin.readline()
				skip -= 1
			
			for szl in fin:
				aln = szl.strip().split("\t")
				self.mapq[aln[0]] = aln[1]
		
		return 
			
			
class MultiRunLog(object):
	def __init__(self):
		self.runlogs = []
		self.stat_order = []
		self.num_logs = 0
	
	def add_log(self, rlog):
		self.runlogs.append(rlog)
		self.num_logs += 1
		
		if len(self.stat_order)==0:
			self.stat_order = list(rlog.stat_order)
	
	def get_stats_table(self):
		
		if self.num_logs==0:
			return None
		
		outlines = []
		# append header
		lout = ["sample"] + self.stat_order
		outlines.append(lout)
		
		# loop through logs
		for rlog in self.runlogs:
			lout = [re.sub("\.bam$", "", rlog.bam)]
			for stat in self.stat_order:
				lout.append(rlog.stats[stat])
			
			outlines.append(map(str, lout))
		
		return outlines
	
	def get_mapq_table(self):
		if self.num_logs == 0:
			return None
			
		# first we have to collect all of the MAPQ values present in all 
		# logs
		mapq = {}
		for rlog in self.runlogs:
			for k in rlog.mapq.keys():
				mapq[k] = []
		
		# loop through again but this time append the values
		for rlog in self.runlogs:
			for k in rlog.mapq.keys():
				mapq[k].append(rlog.mapq[k])
		
		# make a table
		ksort = map(str, sorted(map(int, mapq.keys())))
		
		outlines = []
		# make header
		lout = ["sample"] + ksort
		outlines.append(lout)
		
		for rlog in self.runlogs:
			lout = [re.sub("\.bam$", "", rlog.bam)]
			for q in ksort:
				if q in rlog.mapq:
					lout.append(rlog.mapq[q])
				else:
					lout.append(0)
			
			outlines.append(map(str, lout))
		
		return outlines
			
		
		
				

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Combine read alignment/assignment stats from multiple 'scrna-gcounts-bd.py' runs.")
parser.add_argument('runlogs', type=str, nargs="+", 
	help="Specify 2 or more 'runlog' files from 'scrna-gcounts-bd.py' output")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

