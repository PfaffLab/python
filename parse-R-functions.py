#!/usr/bin/python
#==============================================================================
# parse-R-functions.py
#
# Shawn Driscoll
# 20161104
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Dumps all functions from a single .R file into individual files named 
# by the function name.
#==============================================================================

import sys, argparse, math, re
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
	flines = []
	findex = []
	lnum = 0
	total_lines = 0
	i = 0
	j = 0
	last_index = False

	#
	# read entire file into a list
	fin = open(args.R_file, "r")
	lnum = -1
	for szl in fin:
		flines.append(szl)
		lnum += 1

		if re.search("^#", szl):
			continue

		#
		# check if this line is a function definition line
		if re.search("\<\- function", szl):
			#sys.stderr.write(szl)
			if re.search("^\s", szl):
				# not a root function call
				sys.stderr.write("....NOT a root function definition\n")
			else:
				# root function call
				findex.append([parse_fname(szl), lnum, lnum])

	fin.close()
	
	total_lines = lnum+1	

	# 
	# run back through the index and find the limits of each function definition
	for j in range(len(findex)):
		ll = findex[j]
		
		# first check to see if there are any # lines above it
		i = ll[1]-1
		while i > 0:
			if re.search("^#", flines[i]):
				ll[1] -= 1

			if re.search("^}", flines[i]):
				break

			if re.search("^\s", flines[i]):
				break

			i -= 1

		i = ll[2]+1
		while i < total_lines:
			if re.search("^}", flines[i]):
				ll[2] = i
				break

			i += 1
			if i==total_lines:
				ll[2] = i

		fout_stub = re.sub("\.", "_", ll[0])

		fout_full = "{}/{}.R".format(args.o, fout_stub)
		fout = open(fout_full, "w")

		# write source file as first line
		fout.write("# Shawn Driscoll, sdriscoll@salk.edu\n# Originally in {}\n".format(args.R_file))

		sys.stderr.write("--- new function {} ---\n".format(fout_stub))
		i = ll[1]
		while i <= ll[2]:
			fout.write(flines[i])
			i += 1
		
		fout.close()

	return 0


def parse_fname(sz):
	lparts = sz.split()
	return(lparts[0])

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Dumps all functions from a single .R file into individual files named by the function name.")
parser.add_argument('R_file', type=str, help="R code file.")
parser.add_argument("-o", type=str, default="./", 
	help="Output path to dump individual .R files [current directory]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

