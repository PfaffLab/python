#!/usr/bin/python
#==============================================================================
# prep-bd-for-kallisto.py
#
# Shawn Driscoll
# 20170801
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# After running 'prep-bd-reads' and then demultiplexing the cells run this
# script to pull the UMI barcodes out to a second file so that the reads
# can be processed with 'kallisto pseudo --umi'
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
import subprocess as sp

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
	left = []
	lcount = 0
	rcount = 0
	fout_name = None
	
	if not isfile(args.fastq):
		error_message("Input file {} does not exist".format(args.fastq))
		
		
	tmp = (args.fastq).split(".")
	if re.search("^\.", args.fastq):
		fout_name = ".{}.umi".format(tmp[1])
	else:
		fout_name = "{}.umi".format(tmp[0])
	#
	# read both files and handle each read one at a time
	#

	if re.search("\.gz$", args.fastq):
		# need to gunzip this guy
		p1 = sp.Popen("gunzip -c {}".format(args.fastq).split(), stdout=sp.PIPE)
		fin1 = p1.stdout
	else:
		fin1 = open(args.fastq, "r")

	fout = open(fout_name, "w")
	message("writing results to {}".format(fout_name))
	# loop through the fastq file
	for szl in fin1:
		lcount += 1
		szl = szl.strip()
		left.append(szl)
				
		if lcount==4:
			rcount += 1

			if (rcount % 1e5) == 0:
				progress_message("parsed {} reads".format(rcount))

			# stop and deal with the read. parse out what we need from the 
			# first mate
			tmp = left[0].split(":")
			umi = tmp[-1]
			fout.write(umi)
			fout.write("\n")
						
			lcount = 0
			left = []

	fout.close()
	fin1.close()
	
	progress_message("parsed {} reads".format(rcount))	
	sys.stderr.write("\n")

	return 0

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))

def progress_message(sz):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r[{}] {}".format(time_string(), sz))

def message_mp(sz, name, lock):
	lock.acquire()
	sys.stderr.write("[{}] {}\n".format(name, sz))
	lock.release()
	return


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Prep BD Precise single cell reads by moving the well barcode from the first mate to the end of the second mate read and moving the MI barcode to the read name.")
parser.add_argument('fastq', type=str, help="Post demultiplexed FASTQ with UMI barcode at the end of the read name field")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

