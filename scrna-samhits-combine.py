#!/usr/bin/python
#==============================================================================
# scrna-samhits-combine.py
#
# Shawn Driscoll
# 20180208
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Load output from multiple cells that have been quantified with scrna-samhits.py.
# all that means is the first line is "#total_reads=X" statement and the 
# following rows are tab-delim "gene_name\tumi_count".
# This script produces that matrix output that kalliso pseudo puts out or 
# how cellranger stores the counts on disk. 
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime, time
from Basics import messages as ms

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
COUNTS = "matrix.tsv"
GENES = "matrix.genes"
SAMPLES = "matrix.samples"

#==============================================================================
# main
#==============================================================================

def main(args):

	if not isfile(args.batch_file):
		ms.error_message("Input batch file does not exist ({})".format(args.batch_file))
		return 1
	
	flist = []
	err_flag = False
	
	ms.message("Checking batch file")
	t0 = time()
	
	# load the batch file and confirm all of the input files exist.
	with open(args.batch_file, "r") as fin:
		for szl in fin:
			fname = szl.strip()
			if not isfile(fname):
				err_flag = True
				ms.error_message("{} does not exist".format(fname))
			else:
				flist.append(fname)
	
	ms.time_diff(t0)
	
	if err_flag:
		return 1
	
	rres = core(flist)

	return 1

##
# core processing. loop through files in flist to build target table and 
# hit table
def core(flist):
	
	dtargets = {}
	ltargets = []
	target_idx = 0
	dsamples = {}
	lsamples = []
	sample_idx = 0
	lhits = []
	
	ms.message("Loading hits from {} files".format(len(flist)))
	t0 = time()
	for f in flist:
		sample_name = f.split(".")[0]
		lsamples.append(sample_name)
		# keep track of total reads
		dsamples[sample_name] = 0
		
		with open(f, "r") as fin:
			for szl in fin:
				szl = szl.strip()
				if szl[0]=="#": 
					r = re.search("^\#total\_reads\=([0-9]+)", szl)
					if r:
						dsamples[sample_name] = r.group(1)
					
					continue
				
				aln = szl.split("\t")
				
				if aln[0] not in dtargets:
					dtargets[aln[0]] = target_idx
					ltargets.append(aln[0])
					target_idx += 1
				
				if float(aln[1]) > 0:
					# record the hit as target index, sample index, hit count
					lhits.append([dtargets[aln[0]], sample_idx, aln[1]])
		
		# increment sample index
		sample_idx += 1
	
	ms.time_diff(t0)
		
	##
	## finished parsing files
	##
	
	ms.message("Writing results to your base.")
	t0 = time()
	
	# create count output
	with open(COUNTS, "w") as fout:
		for l in lhits:
			sz = "\t".join(map(str, l))
			fout.write(sz + "\n")
	
	# create samples file
	with open(SAMPLES, "w") as fout:
		for sid in lsamples:
			l = [sid, dsamples[sid]]
			fout.write("\t".join(map(str, l)))
			fout.write("\n")
	
	# create genes file
	with open(GENES, "w") as fout:
		for tid in ltargets:
			fout.write(tid + "\n")
	
	ms.time_diff(t0)
	
	return 0


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('batch_file', type=str, help="Batch input file with one line per UMI count file.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

