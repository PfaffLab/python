#!/usr/bin/env python
#==============================================================================
# make-genome-bin-index.py
#
# Shawn Driscoll
# 20160412
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script uses a 'fai' index for a reference or any tsv format with 
# reference name in the first column and length in the second.
#==============================================================================

import sys, argparse
#import math
#import re
from os.path import isfile, expanduser
#import hashlib
#import subprocess as sp
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
	rname = ""
	left = 0
	right = 0
	w = args.w

	# open index
	fin = open(args.ref_sizes, "r")

	for szl in fin:
		aln = szl.strip().split("\t")
		rname = aln[0]

		left = 0
		right = left + w
		while right < int(aln[1]):
			print "\t".join(map(str, [rname, left, right]))
			left += w
			right = left + w

		# if length is not evenly divisible by the window size then there will 
		# be one last window
		if (float(aln[1]) % w) != 0:
			print "\t".join(map(str, [rname, left, aln[1]]))


	fin.close()


	return 0



#==============================================================================
# defs
#==============================================================================




#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('ref_sizes', type=str, help="Reference names and lengths in tab delim format")
parser.add_argument('-w', type=int, default=25, 
	help="Bin window width [25]")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

