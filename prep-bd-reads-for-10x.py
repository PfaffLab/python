#!/usr/bin/python
#==============================================================================
# prep-bd-reads.py
#
# Shawn Driscoll
# 20170807
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script converts the BD reads to a format compatible with the 
# pachter lab 10x pipeline for splitting into cells, etc.
# in this format each read pair is a barcode file and a read file. the read
# file contains the reads and umis interleaved as READ/UMI/READ/UMI etc
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
import subprocess as sp
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
	left = []
	right = []
	lcount = 0
	rcount = 0

	index_format = "read-I1_si-ATCGCTCC_lane-001-chunk-001.fastq"
	read_format = "read-RA_si-ATCGCTCC_lane-001-chunk-001.fastq"
	
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
	

	with open(index_format, "w") as fout_index, open(read_format, "w") as fout_reads:
			
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

				# parse out the cell barcode and the umi
				
				left_read = left[1]
				left_qc = left[3]

				well_bc = left_read[0:8]
				well_qc = left_qc[0:8]
				
				mi = left_read[8:16]
				mi_qc = left_qc[8:16]

				#
				# write the cellbarcode to the index fastq
				#
				fout_index.write(left[0]+"\n")
				fout_index.write(well_bc+"\n")
				fout_index.write("+\n")
				fout_index.write(well_qc+"\n")

				#
				# write the read to the reads file
				#
				fout_reads.write(right[0]+"\n")
				fout_reads.write(right[1]+"\n")
				fout_reads.write(right[2]+"\n")
				fout_reads.write(right[3]+"\n")

				#
				# write the UMI
				#
				fout_reads.write(right[0]+"\n")
				fout_reads.write(mi+"\n")
				fout_reads.write("+\n")
				fout_reads.write(mi_qc+"\n")

				lcount = 0
				right = []
				left = []

	fin1.close()
	fin2.close()
	
	progress_message("parsed {} reads".format(rcount))	
	sys.stderr.write("\n")

	cmd = "gzip -f {}".format(index_format)
	os.system(cmd)
	cmd = "gzip -f {}".format(read_format)
	os.system(cmd)

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
parser.add_argument('-b', type=int, default=8, action="store", 
	help="Cell barcode length. 14 or 16 for 10x and 8 for BD precise [8]")
parser.add_argument('-u', type=int, default=8, action="store", 
	help="UMI barcode length. 8 for BD Precise and 10 for 10x. [8]"

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

