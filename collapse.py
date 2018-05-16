#!/usr/bin/env python
#==============================================================================
# collapse.py
#
# Shawn Driscoll
# 20130617
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Used to collapse counts by some id column (as specified in the command line) 
#==============================================================================

import sys, argparse

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	hid = {}
	lheader = []
	num_info_cols = args.first_numeric-1

	# check input file
	if not file_exists(args.infile):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.infile)
		return 1


	fin = open(args.infile, "r")

	# deal with the header row.  pop out the collapse column name and prepend it 
	# to the front so that it becomes the first column
	szl = fin.readline()
	lheader = szl.strip().split("\t")
	col_name = lheader.pop(args.collapse_col)
	lheader.insert(0, col_name)


	for szl in fin:
		ll = szl.strip().split("\t")

		tid = ll[args.collapse_col]
		if tid not in hid:
			hid[tid] = {}
			hid[tid]['info'] = [set([]) for i in range(num_info_cols)]
			hid[tid]['values'] = []

		# update info fields
		idx = 0
		for i in range(args.first_numeric):
			if i != args.collapse_col:
				hid[tid]['info'][idx].update([ll[i]])
				idx += 1

		# update values
		if len(hid[tid]['values']) == 0:
			hid[tid]['values'] = list(map(float, ll[args.first_numeric:]))
		else:
			idx = 0
			for i in range(args.first_numeric, len(ll)):
				# sum or max, depending on -m setting
				if args.m:
					# use max value
					hid[tid]['values'][idx] = max(float(ll[i]), hid[tid]['values'][idx])
				else:
					# sum
					hid[tid]['values'][idx] += float(ll[i])
					
				idx += 1


	#---------- print results
	sys.stdout.write("\t".join(lheader) + "\n")
	for tid in sorted(hid.keys()):
		sys.stdout.write(tid + "\t")
		for i in range(num_info_cols):
			ltmp = list(hid[tid]['info'][i])
			sys.stdout.write(",".join(ltmp))
			sys.stdout.write("\t")

		sys.stdout.write("\t".join(map(str, hid[tid]['values'])))
		sys.stdout.write("\n")


	fin.close()

	return 0


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

def is_number(s):
	result = False

	try:
		float(s)
		result = True
	except ValueError:
		result = False

	return result

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Used to collapse counts by some id column (as specified in the command line).")
parser.add_argument('infile', type=str, help="Input file")
parser.add_argument('collapse_col', type=int, help="Column used to collapse file. This column will become the first column in the output file.")
parser.add_argument('first_numeric', type=int, help="First numeric column for summing values in collapsed rows")
parser.add_argument('-m', dest="m", action="store_const", const=True, default=False, help="Report maximum value for collapsed rows (default sum)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
