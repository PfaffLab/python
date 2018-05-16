#!/usr/bin/env python
#==============================================================================
# mutate-genome.py
#
# Shawn Driscoll
# 20130403
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Provided some error rate this script reads in FASTA file and introduces
# single point mutations. A BED format file is also generated containing all of 
# the introduced mutations.
#==============================================================================

import sys, argparse, math
from random import gauss, random, sample
from scipy.stats import norm
import numpy as np
import numpy.random as npr


#==============================================================================
# main
#==============================================================================
def main(args):

	# variables
	seq_header = ""
	seq = ""
	prob = 0
	ppos = 0
	lzl = []
	alpha = 0
	n = 0
	i = 0
	ref_name = ""

	# 
	# check input file
	#
	if args.infile == "stdin" or args.infile == "-":
		fin = sys.stdin
	else:
		if not file_exists(args.infile):
			sys.stderr.write("Error: input file doesn't exist!\n")
			return(1)

		#
		# read in fasta and print it back out as we go
		fin = open(args.infile, "r")
	
	# fasta output
	fout = open(args.stub + ".fa", "w")
	# mutation dictionary output
	fsnps = open(args.stub + ".snp.bed", "w")

	#
	# copy variable for no reason at all
	alpha = args.error_rate

	for szl in fin:

		if szl[0] == ">":
			fout.write(szl)
			# snag the reference name
			ref_name = szl[1:].split(" ")[0] # only in python, bro
			ppos = 0
		else:
			szl = szl.strip()

			# make a probabilty vector for this row
			# x = abs(npr.randn(len(szl)))
			# y = x * -1
			# prob = 1 - (norm.cdf(x) - norm.cdf(y))
			prob = npr.rand(len(szl))
			prob_ci = np.where(prob < alpha)[0]

			if len(prob_ci) > 0:
				n = len(prob_ci)
				i = 0
				lzl = list(szl)
				while i < n:
					ref_base = lzl[prob_ci[i]]
					lzl[prob_ci[i]] = mutate_base(ref_base)
					fsnps.write("{:s}\t{:d}\t{:d}\tref={:s};alt={:s}\n".format(ref_name, prob_ci[i]+ppos, prob_ci[i]+ppos+1, ref_base, lzl[prob_ci[i]]))
					i += 1

				szl = "".join(lzl)

			fout.write(szl + "\n")

			ppos += len(szl)


	fout.close()
	fin.close()
	fsnps.close()

	return 0


#==============================================================================
# defs
#==============================================================================

#
# generate_base_error
#
def mutate_base(base):
	bases_uc = ['A', 'C', 'T', 'G']
	bases_lc = ['a', 'c', 't', 'g']
	base_use = None

	if base == "N" or base == "n":
		return(base)

	e_base = ""

	if base in bases_uc:
		base_use = bases_uc
	else:
		base_use = bases_lc

	e_base = sample(base_use, 1)[0]
	while e_base == base:
		e_base = sample(base_use, 1)[0]

	return(e_base)


#
# check for file
#
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

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Generates a mutated copy of the input FASTA and a report of the mutations as a bed file.")
parser.add_argument('infile',type=str,help="FASTA input file")
parser.add_argument('-o', dest="stub", default="reads", type=str, help="Output reads file stub (default: reads)")
parser.add_argument("-e", dest="error_rate", default=0, type=float, help="Base error rate. Illumina typical could be 0.001 (default: 0)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
