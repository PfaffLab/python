#!/usr/bin/python
#==============================================================================
# kmer-transcript-library.py
#
# Shawn Driscoll
# 20170315
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Extracts each transcript from a FASTA and generates a single kmer
# file per feature. 
#==============================================================================

import sys, argparse
# import math, re
from os.path import isfile, expanduser, isdir
from os import system

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
	szfai = "{}.fai".format(args.fasta)
	rres = None
	tid2gid = {}
	tid2gname = {}
	has_info = False
	fout_name = "{}.tmapp".format(args.fasta)
	ntid = 0

	#
	# check for fasta index
	if not fai_exists(args.fasta):
		# need to build on
		message("building FAI index")
		build_fai(args.fasta)

	#
	# parse the fai
	dfai = parse_fai(szfai)

	#
	# check for output folder
	if not isdir(args.o):
		message("Creating output folder {}".format(args.o))
		runcmd("mkdir {}".format(args.o))

	#
	# test this idea with the first few
	for tid in dfai.keys()[0:10]:
		#
		# fetch sequence and kmer it
		message("Processing {}".format(tid))
		cmd = "samtools faidx {} {} > tmp.fa".format(args.fasta, tid)
		runcmd(cmd, verbose=False)
		cmd = "kmercountexact.sh in=tmp.fa k={} rcomp=f out={}/{}.k{}.fa overwrite=t 2>kmer.log".format(args.kmer_length, args.o, tid, args.kmer_length)
		runcmd(cmd, verbose=False)
		ntid += 1

	# print results
	#rres = print_results(dtid, fout_name)

	return 0


#==============================================================================
# defs
#==============================================================================

def parse_fai(fai):
	dfai = {}
	fin = open(fai, "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		dfai[aln[0]] = int(aln[1])

	fin.close()
	return dfai

#
# check for fai index for the supplied fasta
def fai_exists(fa):
	return(isfile("{}.fai".format(fa)))

#
# if the supplied fasta file doesn't have an fai index then this function
# calls samtools to build it
def build_fai(fa):
	cmd = "samtools faidx {}".format(fa)
	runcmd(cmd)
	return 0

#
# run a system level command. used for running alignment and samtools commands
def runcmd(cmd, verbose=True):
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	rres = system(cmd)
	return rres

#
# print a message to STDERR. appends the newline for you.
def message(sz):
	sys.stderr.write("[kmer-transcript-library] " + sz + "\n")


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")

# -- required args
parser.add_argument('fasta', type=str, help="Input FASTA")

# -- options
parser.add_argument('-o', default="kmerLib", type=str, action="store", 
	help="Output folder for kmer library [kmerLib]")
parser.add_argument('-k', '--kmer-length', default=31, action="store", 
	help="Kmer length for matching. Must be an odd number [31]")

# -- parse command line
args = parser.parse_args()


if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

