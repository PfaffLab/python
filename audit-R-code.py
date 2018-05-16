#!/usr/bin/python
#==============================================================================
# audit-R-code.py
#
# Shawn Driscoll
# 20161104
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Finds all functions that are called within an R file.
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
	findex = {}
	fcalls = {}
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

			fname = parse_fname(szl.strip())
			if fname not in findex:
				findex[fname] = []

			findex[fname].append(lnum)


			#sys.stderr.write(szl)
#			if re.search("^\s", szl):
#				# not a root function call
#				sys.stderr.write("....NOT a root function definition\n")
#			else:
#				# root function call
#				fname = parse_fname(szl)
#				if fname not in findex:
#					findex[fname] = lnum
				

	fin.close()
	
	# loop through all lines and find function calls
	for szl in flines:
		rres = re.findall("([a-zA-Z0-9\_\.]+)\(", szl)
		if len(rres) > 0:
			for term in rres:
				if term not in fcalls:
					fcalls[term] = 0
				fcalls[term] += 1


	for term in sorted(fcalls.keys()):
		lout = [term, fcalls[term], "."]

		if term in findex:
			# function is defined in this file:
			lout[2] = "defined @ {}".format(",".join(map(str,findex[term])))

		elif args.p is not None:
			#
			# check path for files that match function call names. print out only those 
			# functions.

			fname = "{}/{}.R".format(args.p, reformat_fname(term))
			if isfile(fname):
				# got it
				lout[2] = fname

		print "\t".join(map(str, lout))

	return 0


def parse_fname(sz):
	lparts = sz.split()
	return(lparts[0])

def reformat_fname(sz):
	sz_hat = re.sub("\.", "_", sz)
	return sz_hat

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Scans an R code file for function calls. Checks if the functions are defined within the file or if they are missing.")
parser.add_argument('R_file', type=str, help="R code file.")
parser.add_argument("-p", type=str, default=None, 
	help="Path to folder of R source files for functions dumped out to one funciton per file (from parse-R-functions.py)")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

