#!/usr/bin/python
#==============================================================================
# smg-parse.py
#
# Shawn Driscoll
# Jan 24, 2017
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse the output of 'seal-mappability-graph.py' to create a per-transcript
# report of mappability and which other transcripts it shares sequence with.
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
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
	dmap = {}
	dinfo = {}

	# 
	# open fai
	sys.stderr.write("Loading {}\n".format(args.fasta_fai))
	fin = open(args.fasta_fai, "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		dmap[aln[0]] = []

	fin.close()

	#
	# open input file and load a dict up with the 'reference' and any 
	# 
	sys.stderr.write("Loading {}\n".format(args.smg_output))
	fin = open(args.smg_output, "r")

	# skip the header
	szl = fin.readline()
	for szl in fin:
		aln = szl.strip().split("\t")

		if float(aln[2]) > args.t:
			dmap[aln[0]].append([aln[1], float(aln[2])])

		if aln[0] not in dinfo:
			dinfo[aln[0]] = [float(aln[3]), int(aln[4]), float(aln[5])]

	fin.close()

	#
	# analyze the groups
	#
	sys.stderr.write("Processing graph and writing report\n")
	for tid in sorted(dmap.keys()):
		num_hits = len(dmap[tid])

		lout = [tid, num_hits, 0, "u", "u"]

		if num_hits > 1:

			# sort the hits
			ssim = [x[1] for x in dmap[tid]]
			ssim_order = list(np.argsort(ssim))

			# highest sim first
			lout[2] = dmap[tid][ssim_order[-1]][1]
			lout[3] = dmap[tid][ssim_order[-1]][0]


		elif num_hits > 0:
			# only 1
			lout[2] = dmap[tid][0][1]
			lout[3] = dmap[tid][0][0]

		print "\t".join(map(str, lout))

	sys.stderr.write("done\n")

	return 0



#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Parse the output of 'seal-mappability-graph.py' to create a per-transcript report of mappability and which other transcripts it shares sequence with.")
parser.add_argument('fasta_fai', type=str, help="FASTA FAI for the reference used in 'seal-mappability-graph.py'. Used to establish the base set of all transcripts.")
parser.add_argument('smg_output', type=str, help="Input file")
parser.add_argument('-t', type=float, default=0.01, 
	help="Threshold for similarity. Only connections above this threshold are retained. [0.05]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

