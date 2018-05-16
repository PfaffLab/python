#!/usr/bin/python
#==============================================================================
# asi-miso-index-from-indexG.py
#
# Shawn Driscoll
# 20170503
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Builds a MISO index from the MISO GFF type output of 'asi-build-indexG'
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
import os

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
	aln = []
	szl = ""
	fname = ""
	fin = None
	htypes = {}
	k = ""
	childs = []
	pid = None
	tmp = None
	nthreads = 2
	nrunning = 0
	
	#
	# open gff and split it into sub-gffs
	#
	fin = open(args.gff, "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		k = aln[1]
		if k not in htypes:
			htypes[k] = []
		
		# push the row from the file into the dict
		htypes[k].append(szl)
	
	fin.close()
	
	#
	# write new gff files for each type
	#
	for k in sorted(htypes.keys()):
		# create output file
		fname = "{}_{}".format(k, args.gff)
		fout = open(fname, "w")
		# write lines to file
		for szl in htypes[k]:
			fout.write(szl)
		# close it
		fout.close()
		
		#
		# build index
		#
		pid = os.fork()
		if pid==0:
			run_miso_index(fname, k)
			os._exit(0)
		else:
			childs.append(pid)
			nrunning += 1
		
		if nrunning >= nthreads:
			# wait
			for tmp in childs:
				os.waitpid(tmp, 0)
			
			childs = []
			nrunning = 0
		
	if nrunning > 0:
		for tmp in childs:
			os.waitpid(tmp, 0)

	return 0

#==============================================================================
# defs
#==============================================================================

def run_miso_index(fname, outdir):
	cmd = "index_gff --index {} {} 2>/dev/null 1>/dev/null".format(fname, outdir)
	return runcmd(cmd)

#
# run a system level command. used for running alignment and samtools commands
def runcmd(cmd, verbose=True):
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	rres = os.system(cmd)
	return rres

#
# print message to stderr
def message(sz):
	sys.stderr.write("[MAIN] " + sz + "\n")


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Builds a MISO index from the MISO GFF type output of 'asi-build-indexG'")
parser.add_argument('gff', type=str, help="GFF3 output from 'asi-build-indexG'")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

