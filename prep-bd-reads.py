#!/usr/bin/python
#==============================================================================
# prep-bd-reads.py
#
# Shawn Driscoll
# 20170615
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script takes paired end reads from a BD precise scRNA run and 
# generates a new single-end file containing the second-mate reads 
# annotated with the well and MI barcodes from the first mate. both barcodes
# are added to the read name (MI last) and the well barcode is added to the
# read so that we can demultiplex the wells with fastx_barcode_splitter. 
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
	right = []
	lcount = 0
	rcount = 0
	fout_name = "{}.fastq".format(args.o)
	
	if not isfile(args.mate1):
		error_message("Input file {} does not exist".format(args.mate1))
	
	if not isfile(args.mate2):
		error_message("Input file {} does not exist".format(args.mate2))
	
	#
	# read both files and handle each read one at a time
	#

	if re.search("\.gz$", args.mate1):
		# need to gunzip this guy
		p1 = sp.Popen("gunzip -c {}".format(args.mate1).split(), stdout=sp.PIPE)
		fin1 = p1.stdout
	else:
		fin1 = open(args.mate1, "r")
	
	if re.search("\.gz$", args.mate2):
		# need to gunzip this guy
		p2 = sp.Popen("gunzip -c {}".format(args.mate2).split(), stdout=sp.PIPE)
		fin2 = p2.stdout
	else:
		fin2 = open(args.mate2, "r")
	

	fout = open(fout_name, "w")
		
	# make loop runs over first mate file
	for szl in fin1:
		lcount += 1
		szl = szl.strip()
		left.append(szl)
		szl2 = (fin2.readline()).strip()
		right.append(szl2)
				
		if lcount==4:
			rcount += 1

			if (rcount % 1e5) == 0:
				progress_message("parsed {} reads".format(rcount))

			# stop and deal with the read. parse out what we need from the 
			# first mate
			left_read = left[1]
			left_qc = left[3]
			well_bc = left_read[0:8]
			well_qc = left_qc[0:8]
			mi = left_read[8:16]
			
			#
			# setup the second mate read with barcodes
			#
			
			# split off anything after an initial whitespace character
			tmp = right[0].split()
			# append both the well barcode and the MI barcode to the read name
			rname = tmp[0] + ":{}:{}".format(well_bc, mi)
			# append the well barcode to the read and the barcode quals to the quals of the read
			right[1] += well_bc
			right[3] += well_qc
			# set the name
			right[0] = rname
			
			# print augmented second mate out to the new file
			for l in right:
				fout.write(l + "\n")
			
			lcount = 0
			right = []
			left = []

	fout.close()
	fin1.close()
	fin2.close()
	
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
parser.add_argument('mate1', type=str, help="First mate file with barcodes")
parser.add_argument('mate2', type=str, help="Second mate file.")
parser.add_argument('-o', type=str, default="bdPrepped", 
	help="Stub for output FASTQ file")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

