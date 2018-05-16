#!/usr/bin/env python
#==============================================================================
# seal-mappability-reduce.py
#
# Shawn Driscoll
# June 2016
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script is meant to run on the output of seal-mappability-graph.py. 
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib

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

#==============================================================================
# main
#==============================================================================

def main(args):

	hpares = {}
	lnomatch = []

	fin = open(args.file, "r")
	# skip header
	tmp = fin.readline()

	for szl in fin:
		aln = szl.strip().split("\t")

		if aln[1] == "none":
			lnomatch.append(szl)
			continue

		k = make_keys(aln)

		if k[0] in hpares or k[1] in hpares:
			if k[0] in hpares:
				hpares[k[0]].append(float(aln[2]))
			elif k[1] in hpares:
				hpares[k[1]].append(float(aln[2]))
		else:
			hpares[k[0]] = [float(aln[2])]

	fin.close()

	print "tid1\ttid2\tsim1\tsim2\thmean"

	for k in hpares.keys():
		t = k.split("|")
		tid = t[0:2]
		print sorted(tid)
		print hpares[k]
		hmean = harmonic_mean(hpares[k])
		print "\t".join([t[0], t[1], str(hpares[k][0]), str(hpares[k][1]), 
			str(hmean)])	
#		print "\t".join([t[0], t[1], str(hmean)])
#		print "\t".join([t[1], t[0], str(hmean)])

	return 0

def make_keys(v):
	k1 = "{}|{}".format(v[0], v[1])
	k2 = "{}|{}".format(v[1], v[0])

	return(list([k1, k2]))

def harmonic_mean(v):

	n = len(v)
	vsum = 0
	vmean = 0

	for i in range(n):
		vsum += 1.0/v[i]

	vmean = vsum/float(n)

	return(1.0/vmean)

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Reduce the mappability table down to unique pairs")
parser.add_argument('file', type=str, help="Output of seal-mappability-graph.py")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

