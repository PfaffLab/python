#!/usr/bin/env python
#==============================================================================
# template.py
#
# Shawn Driscoll
# 20151207
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Binning mpileup into fixed windows.
#==============================================================================

import sys, argparse, re
from os.path import isfile
import subprocess as sp

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

	if not isfile(args.bam):
		sys.stderr.write("Error: input BAM doesn't exist ({})".format(args.bam))
		return(1)

	# check args
	if args.r is not None:
		if args.f is None:
			sys.stderr.write("If specifying a region then you have to also use -f\n")
			return(1)

	try:
		href_lengths = get_ref_sizes(args.bam)
		href_num_bins = num_bins(href_lengths, args.b)
		href_bins = prepare_bin_vectors(href_num_bins)
	except Exception as e:
		raise

	try:
		href_bins_counted = count_bins(args.bam, href_bins, args.b, args.q, args.f, args.r)
	except Exception as e:
		raise

	print_results(href_bins_counted, args.b)

	return 0

def print_results(countedBins, binSize):
	# print results out

	for k in sorted(countedBins.keys()):
		for i in range(len(countedBins[k])):
			print "{}\t{}\t{}\t{:f}".format(k, i*binSize, countedBins[k][i], countedBins[k][i]/float(binSize))

def count_bins(bam, hbinvec, binSize, mapq, fasta, region):
	# count per-base abundance per bin
	binSize = int(binSize)
	szl = ""
	aln = []
	bin = 0

	cmd = "samtools mpileup -q {}".format(mapq)

	if fasta is not None:
		cmd += " -f {}".format(fasta)

		if region is not None:
			cmd += " -r {}".format(region)

	cmd += " {}".format(bam)

	sys.stderr.write("CMD: {}\n".format(cmd))

	# call process
	p1 = sp.Popen(cmd.split(), stdout=sp.PIPE)
	for szl in p1.stdout:
		aln = szl.strip().split("\t")
		bin = int(aln[1])/binSize
		hbinvec[aln[0]][bin] += float(aln[3])

	return(hbinvec)


def get_ref_sizes(bam):
	# open the bam and get the header. parse that out to collect
	# the lengths of each reference

	p1 = sp.Popen("samtools view -H {}".format(bam).split(), stdout=sp.PIPE)
	hrefs = {}
	aln = []
	atmp = []
	refid = ""

	for szl in p1.stdout:
		szl = szl.strip()
		if re.search("^@SQ", szl):
			aln = szl.split("\t")
			atmp = aln[1].split(":")
			refid = atmp[1]
			atmp = aln[2].split(":")
			hrefs[refid] = int(atmp[1])

	return(hrefs)

def num_bins(hlens, binSize):
	# figure out total bins per reference 
	binSize = int(binSize)
	hnumBins = {}

	for k in hlens.keys():
		hnumBins[k] = hlens[k]/binSize + 1

	return(hnumBins)

def prepare_bin_vectors(hbins):
	# make bin vectors for each reference based on the number of bins dict

	hbin_vectors = {}
	for k in hbins.keys():
		hbin_vectors[k] = [0 for i in range(hbins[k])]

	return(hbin_vectors)

def call_mpileup(bam):
	return 0

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Runs samtools mpileup and bins the results in fixed windows by summing total bases. Results are printed to stdout")
parser.add_argument('bam', type=str, help="Coordinate sorted alignments")
parser.add_argument('-r', type=str, action="store", default=None, 
	help="Region to return pileup (see samtools mpileup help)")
parser.add_argument('-f', type=str, action="store", default=None, 
	help="Path to fai indexed reference fasta (required if using -r)")
parser.add_argument('-q', type=int, action="store", default=1, 
	help="Minimum MAPQ to count [1]")
parser.add_argument('-b', type=int, action="store", default=16000, 
	help="Bin size in bp [16000]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

