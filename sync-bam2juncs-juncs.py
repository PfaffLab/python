#!/usr/bin/python
#==============================================================================
# sync-bam2juncs-juncs.py
#
# Shawn Driscoll
# 20160819
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This is used to synchronize several bam2juncs junction outputs so they all 
# have the same sorting and the same set of junctions
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

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
	d_jset = {}
	num_files = len(args.flist)

	# must have at least two files
	if num_files < 2:
		sys.stderr.write("Error: you should provide at least two files")
		return 1

	# first pass is to get all junctions
	

	return 0



#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('bam2junc_files', type=str, metavar="flist", nargs="+", 
	help="bam2juncs input files to synchronize")
parser.add_argument("-s", type=str, default="sync", action="store", 
	help="Suffix to append to the input file name [sync]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

