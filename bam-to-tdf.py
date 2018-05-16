#!/usr/bin/python
#==============================================================================
# bam-to-tdf.py
#
# Shawn Driscoll
# 20151204
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script wraps up the process of making a TDF file from 
# a sorted BAM. 
#==============================================================================

import sys, argparse
import subprocess as sp
from os.path import expanduser
from os.path import isfile

# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

_HOME = expanduser("~")
_TDF_UTIL = "/opt/igv/IGVTools/igvtools count"

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	# make sure the bam file exists
	if not isfile(args.bam):
		sys.stderr.write("Input BAM file doesn't exist {}\n".format(args.bam))
		return(1)

	if not isfile(args.genome):
		sys.stderr.write("Genome file does not exist ({})\n".format(args.genome))

	# make sure the bam index is present
	try:
		confirm_bam_index(args.bam)
	except:
		sys.stderr.write("Failed to check for BAM index\n")
		return(1)

	# create the output file name
	s1 = (args.bam).split("/")
	s2 = s1[-1].split(".")
	outputfile = s2[0] + ".tdf"

	# ok build the command
	cmd = _HOME + _TDF_UTIL + " -w {} -z {} -f {} -e {}".format(args.w, args.z, args.f, args.e)
	cmd = cmd + " {} {} {}".format(args.bam, outputfile, args.genome)

	p1 = runcmd(cmd)
	p1.wait()

	return 0

def runcmd(cmd):
	# use Popen to run a system command with no STDOUT or STDERR 
	# capture. return the Popen object
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())
	return p1

def confirm_bam_index(f):
	#
	# check if the bam index file exists

	if isfile("{}.bai".format(f)):
		return 0
	else:
		sys.stderr.write("BAM index is not present - generating one now...\n")
		cmd = "samtools index {}".format(f)
		p1 = runcmd(cmd)
		p1.wait()

	return 0


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Runs IGVTools 'count' command to generate a TDF file an input BAM file.")
parser.add_argument('genome', type=str, help="Genome FASTA")
parser.add_argument('bam', type=str, help="Coordinate sorted BAM alignments relative to the supplied genome")

# igvtools count options
parser.add_argument('-z', type=float, default=7, action="store", 
	help="Maximum zoom level to compute [7]")
parser.add_argument('-w', type=int, default=25, action="store",
	help="Window size over which coverage is averaged [25]")
parser.add_argument('-e', type=int, default=0, action="store", 
	help="The read or feature is extended by the specified distance in bp prior to counting.")
parser.add_argument('-f', type=str, default="mean", action="store", 
	help="A comma separated list specifying window functions to use when reducing the data. Possible values are: min, max, mean, median, p2, p10, p90 and p98. [mean]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

