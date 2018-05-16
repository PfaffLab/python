#!/usr/bin/env python
#==============================================================================
# exactSnp-to-vcf.py
#
# Shawn Driscoll
# 20130618
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Converts exactSnp output to VCF format (or at least a weak VCF format)
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

	# check input file
	if not file_exists(args.infile):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.infile)
		return 1

	fin = open(args.infile, "r")
	szl = fin.readline()

	print "##fileformat=VCFv4.1"
	print "##exactSnp (subread package)"
	print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">"
	print "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">"
	print "##FORMAT=<ID=NL,Number=1,Type=String,Description=\"Noise levels\">"
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.infile

	for szl in fin:
		ll = szl.strip().split("\t")
		if len(ll) == 8:
			if len(ll[3]) > 1:
				# position has multiple alternative alleals
				ltmp = ll[4].split("/")
				alt_hits = 0
				for i in range(len(ltmp)):
					alt_hits += int(ltmp[i])
			else:
				alt_hits = int(ll[4])

			ref_hits = int(ll[5])-alt_hits
			

			lout = [ll[0], ll[1], ".", ll[2], ll[3], ll[6], "."]
			lout.append("DP=%d;DP4=%d,0,%d,0" % (int(ll[5]), ref_hits, alt_hits))
			lout.append("NL")
			lout.append(ll[7])

			print "\t".join(lout)

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
