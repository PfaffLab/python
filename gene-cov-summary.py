#!/usr/bin/env python
#==============================================================================
# gene-cov-summary.py
#
# Shawn Driscoll
# 20130612
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses the output of 'nexpr-coverage' to generate a gene level summary
# of best potential coverage for use in filtering gene lists in 
# downstream analysis.
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
# globals
#==============================================================================

_GID_COL = 1
_COV_COL = 8

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	result = 0

	# check input file
	if not file_exists(args.infile):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.infile)
		return 1

	
	result = core(args)

	return result

def core(args):

	# variables
	szl = ""
	ll = []
	g_hash = {}
	gid = ""
	cov = 0
	stats = [0 for i in range(11)]
	stats_c = list(stats)
	total = 0

	# open input file
	fin = open(args.infile, "r")
	# skip header
	szl = fin.readline()

	for szl in fin:
		ll = szl.strip().split("\t")

		gid = ll[_GID_COL]
		cov = float(ll[_COV_COL])

		if gid not in g_hash:
			g_hash[gid] = 0

		if cov > g_hash[gid]:
			g_hash[gid] = cov

	fin.close()

	for gid in sorted(g_hash.keys()):
		cov = g_hash[gid]
		sys.stdout.write("%s\t%f\n" % (gid, g_hash[gid]))

		if cov == 0:
			stats[0] += 1
		elif cov > 0 and cov <= 0.1:
			stats[1] += 1
		elif cov > 0.1 and cov <= 0.2:
			stats[2] += 1
		elif cov > 0.2 and cov <= 0.3:
			stats[3] += 1
		elif cov > 0.3 and cov <= 0.4:
			stats[4] += 1
		elif cov > 0.4 and cov <= 0.5:
			stats[5] += 1
		elif cov > 0.5 and cov <= 0.6:
			stats[6] += 1
		elif cov > 0.6 and cov <= 0.7:
			stats[7] += 1
		elif cov > 0.7 and cov <= 0.8:
			stats[8] += 1
		elif cov > 0.8 and cov <= 0.9:
			stats[9] += 1
		elif cov > 0.9 and cov <= 1.0:
			stats[10] += 1

		total += 1

	for i in range(len(stats)):
		stats[i] = stats[i]*1.0/total
		if i > 0:
			stats_c[i] = stats[i]+stats_c[i-1]
		else:
			stats_c[i] = stats[i]

	# print stats out to stderr
	sys.stderr.write("Statistics:\n")
	sys.stderr.write("cov\tratio\tcumulative\n")
	sys.stderr.write("0\t%f\t%f\n" % (stats[0], stats_c[0]))
	sys.stderr.write("0 to 0.1\t%f\t%f\n" % (stats[1], stats_c[1]))
	sys.stderr.write("0.1 to 0.2\t%f\t%f\n" % (stats[2], stats_c[2]))
	sys.stderr.write("0.2 to 0.3\t%f\t%f\n" % (stats[3], stats_c[3]))
	sys.stderr.write("0.3 to 0.4\t%f\t%f\n" % (stats[4], stats_c[4]))
	sys.stderr.write("0.4 to 0.5\t%f\t%f\n" % (stats[5], stats_c[5]))
	sys.stderr.write("0.5 to 0.6\t%f\t%f\n" % (stats[6], stats_c[6]))
	sys.stderr.write("0.6 to 0.7\t%f\t%f\n" % (stats[7], stats_c[7]))
	sys.stderr.write("0.7 to 0.8\t%f\t%f\n" % (stats[8], stats_c[8]))
	sys.stderr.write("0.8 to 0.9\t%f\t%f\n" % (stats[9], stats_c[9]))
	sys.stderr.write("0.9 to 1.0\t%f\t%f\n" % (stats[10], stats_c[10]))
	return 0


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('infile', type=str, help="Input file")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
