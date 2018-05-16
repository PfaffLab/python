#!/usr/bin/env python
#==============================================================================
# macs-wrapper.py
#
# Shawn Driscoll
# 20151221
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse a simple text file with names of signal BAMS and provide a 
# input [optionally] to call peaks against. Input can be several files
# merged with samtools [merge] 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
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
# globals
#==============================================================================

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	flist = []

	# open the samples file list and load them into a list
	fin = open(args.flist, "r")
	for szl in fin:
		flist.append(szl.strip())

	fin.close()
	sys.stderr.write("Parsed {} samples from {}\n".format(len(flist), args.flist))
	
	# loop through these and get it on
	for f in flist:
		sid = path_to_sample_id(f)
		# run it
		runmacs(f, args.c, sid)


	return 0


def path_to_sample_id(fpath):
	# split up by /, last part is the file, drop the extension...boom
	
	atmp = fpath.split("/")
	atmp2 = atmp[-1].split(".")
	
	return(atmp2[0])


def runmacs(fsignal, finput, sample):
	
	cmd = "macs14 -t {} -f BAM -g mm --nomodel --shiftsize=150 -n {}".format(fsignal, sample)
	if len(finput) > 0:
		cmd += " -c {}".format(finput)
	sys.stderr.write("CMD: {}".format(cmd))
	pro = sp.Popen(cmd.split())
	pro.wait()
	
	return(0)


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('flist', type=str, help="Text file with one column. Each row should be full path to a signal BAM.")
parser.add_argument('-c', type=str, default="", help="Input BAM for calling peaks against a background")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
#	except Exception as e:
#		sys.stderr.write("\nerror: {}\n".format(e.strerror))

