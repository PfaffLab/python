#!/usr/bin/env python
#==============================================================================
# mirna-parse-hits.py
#
# Shawn Driscoll
# 20160519
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Adapted from juncdb-parse-hits.py. this is used for counting hits of reads
# to a kmer index for micro-rna
#==============================================================================

import sys, argparse 
import math 
import re
import copy
from os.path import isfile, expanduser
import subprocess as sp
import os

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
	dhits = {}

	# parse the hits
	sys.stderr.write("Parsing hits from {}\n".format(args.fa))
	dhits = parse_alignments(args.fa)

	sys.stderr.write("Building final table\n")

	print "\t".join(["name", "unique_hits", "multi_hits", "hits"])

	for tid in sorted(dhits.keys()):

		lout = [tid, dhits[tid]['unique'], dhits[tid]['multi'], dhits[tid]['scaled']]

		print "\t".join(map(str, lout))

	return 0


#==============================================================================
# functions
#==============================================================================

#
# parse_alignments
# this function parses the Seal hits and generates a dict of features hit with 
# counts
def parse_alignments(fname):

	szl = ""
	ltargets = []
	lhits = []
	feature_hits = {}
	wfactor = 1
	num_tied = 0
	idx = 0

	# parse the hits
	fin = open(args.fa, "r")
	for szl in fin:
		if not re.match("^>", szl):
			continue
		
		idx += 1
		if (idx % 1000000) == 0:
			sys.stderr.write("parsed {} aligned reads\n".format(idx))

		wfactor = 1.0
		ltargets, lhits = parse_hit_name(szl)
		
		num_tied = 1
		if len(ltargets) > 1:

			# take best only. hit counts are sorted in descending order already
			hmax = lhits[0]

			# increment through the hit counts to find how many share the sme
			# max match count
			i = 1
			while i < len(lhits):
				if lhits[i] < hmax:
					break

				i += 1

			# hit adjustment factor is 1/<number of best hit features>
			wfactor = 1.0/(i**2)
			num_tied = i

		# count the hits into the features hit adding their names to the dict as we go
		for j in range(0,num_tied):

			if ltargets[j] not in feature_hits:
				feature_hits[ltargets[j]] = { 'unique': 0, 'multi': 0, 'scaled': 0 }

			feature_hits[ltargets[j]]['scaled'] += wfactor
			if num_tied > 1:
				feature_hits[ltargets[j]]['multi'] += 1
			else:
				feature_hits[ltargets[j]]['unique'] += 1

	# final message
	sys.stderr.write("parsed {} aligned reads (done)\n".format(idx))
	
	fin.close()	

	return feature_hits


# ---
# parse a hit name from one of the FASTA lines output from Seal
def parse_hit_name(sz):

	aln = sz.strip().split("\t")

	dhits = {}
	ltargets = []
	lhits = []

	# just one feature means there's only a read name and no targets
	if len(aln)==1:
		return ltargets, lhits

	# loop through targets. because of how the index is made we'll have hits to many mers
	# of the same target so we need to bin them
	for i in range(1,len(aln)):
		alnx = aln[i].split("=")
		alnxName = alnx[0].split("|")

		for i in range(len(alnxName)):
			if alnxName[i] not in dhits:
				dhits[alnxName[i]] = 0

			dhits[alnxName[i]] += int(alnx[1])
	
	# pass to lists
	for tid in dhits.keys():
		ltargets.append(tid)
		lhits.append(dhits[tid])

	# sort by hit count
	if len(lhits) > 1:
		o = list(np.argsort(lhits))
		o.reverse()
		ltmpa = list(lhits)
		ltmpb = list(ltargets)
		for i in range(len(lhits)):
			lhits[i] = ltmpa[o[i]]
			ltargets[i] = ltmpb[o[i]]

	return ltargets, lhits



# --
# runcmd
# run a system level command in subprocess. optionally you can return the process.
# if the process isn't returned then the function waits for the process to finish
def runcmd(cmd, returnProcess=False):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Build junction database for mapping.")
parser.add_argument('fa', type=str, help="seal hits as FASTA")

args = parser.parse_args()

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

