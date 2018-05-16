#!/usr/bin/python
#==============================================================================
# collapse-annotation.py
#
# Shawn Driscoll
# 20180122
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script collapses an annotation table by some column reducing the
# annotation down to one row per unique factor in that column. Values 
# in other columns are collapsed and joined with commas. 
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
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

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	dkey_index = defaultdict(list)
	ltable = []
	lheader = []
	idx = 0
	col_idx = args.c - 1
	num_col = []
	
	if not isfile(args.source_file):
		ms.error_message("Input file does not exist")
		return 1
	
	##
	## open input file and parse it in
	with open(args.source_file, "r") as fin:
		
		if args.H:
			szl = fin.readline()
			lheader = szl.strip().split("\t")
		
		for szl in fin:
			# put each row of the file into a list. also keep track of which rows belong to 
			# which keys in a dict
			aln = szl.strip().split("\t")
			ltable.append(aln)
			num_col.append(len(aln))
			dkey_index[aln[col_idx]].append(idx)
			idx += 1

	##
	# check column count for consistency
	if min(num_col) != max(num_col):
		ms.error_message("Column count is not consistent throughout the input file!")
		return 1

	num_col = min(num_col)

	##
	## build collapsed version
	lout = []
	for kid in sorted(dkey_index.keys()):
		tmp = []
		num_row = len(dkey_index[kid])
		
		for i in range(num_col):
			if i == col_idx:
				continue
			
			if num_row > 1:
				# this key has multiple rows so we have to collapse the value 
				# from the current column within each row into a single string
				ltmp = []
				for j in dkey_index[kid]:
					ltmp.append(ltable[j][i])
				
				# check to see if this field has more than a single level
				set_tmp = set(ltmp)
				if len(set_tmp) == 1:
					# single level
					tmp.append(ltmp[0])
				else:
					# more than one level so keep all of them
					tmp.append(",".join(ltmp))
			
			else:
				# no collapse, single row
				tmp.append(ltable[dkey_index[kid][0]][i])
		
		lout.append("\t".join([kid] + tmp) +  "\n")
	
	if args.H:
		hout = [lheader[col_idx]]
		for i in range(num_col):
			if i != col_idx:
				hout.append(lheader[i])
	
		sys.stdout.write("\t".join(hout))
		sys.stdout.write("\n")
	
	for i in range(len(lout)):
		sys.stdout.write(lout[i]) 

	return 0



#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('source_file', type=str, help="Tab-delim table to collapse")
parser.add_argument('-c', type=int, default=1, 
	help="Column to collapse on (one-based)")
parser.add_argument('-H', action="store_const", const=True, default=False, 
	help="Input file has a header row")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

