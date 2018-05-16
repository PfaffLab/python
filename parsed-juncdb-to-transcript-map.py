#!/usr/bin/env python
#==============================================================================
# parsed-juncdb-to-transcript-map.py
#
# Shawn Driscoll
# 20151211
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Takes output of parsed-junc-db-hits.py and build a transcript centric
# table with each transcript's junction counts, # total, # used
#==============================================================================

import sys, argparse, math, re
from os.path import isfile

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
	tinfo = {}
	tcounts = {}
	tid = ""
	num_j = 0
	used_j = 0
	jcounts = ""


	if not isfile(args.file):
		return 1

	# load the file
	tinfo, tcounts = parse_file(args.file)

	# build output

	# header
	lout = ["transcript_id", "gene_id", "gene_name", "strand", "junction_counts", 
		"num_juncs", "used_juncs"]
	print "\t".join(lout)

	for tid in sorted(tinfo.keys()):
		jcounts = ",".join(map(str, tcounts[tid]))
		num_j = len(tcounts[tid])
		used_j = 0
		for i in range(num_j):
			if tcounts[tid][i] > 0:
				used_j += 1

		lout = tinfo[tid] + [jcounts, num_j, used_j]
		print "\t".join(map(str, lout))

	# done!

	return 0

def parse_file(f):

	# variables
	dtid = {}  # dict of all transcripts
	dtinfo = {}
	dtstrand = {}
	szl = ""
	aln = []
	fin = None
	tid = ""
	ltid = []

	# open file
	fin = open(f, "r")

	# skip header
	szl = fin.readline()
	for szl in fin:
		szl = szl.strip()
		# 3 and 16 are what we want

		aln = szl.split("\t")

		ltid = aln[16].split(",")
		for tid in ltid:
			if tid not in dtid:
				dtinfo[tid] = [tid, aln[15], aln[14], aln[5]]
				dtid[tid] = [] # used for junction counts

			# append so the counts are transcript oriented
			if aln[5] == "+":
				dtid[tid].append(float(aln[3]))
			else:
				dtid[tid] = list([float(aln[3])] + dtid[tid])

	fin.close()

	return (dtinfo, dtid)




#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Takes output of parsed-junc-db-hits.py and build a transcript centric # table with each transcript's junction counts, # total, # used.")
parser.add_argument('file', type=str, help="Parsed juncdb hits from parse-junc-db-hits.py")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

