#!/usr/bin/python
#==============================================================================
# fastq2bam.py
#
# Shawn Driscoll
# 20180226
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Convert FASTQ to BAM. Useful as a way to sort reads by name since 
# samtools is suited for that. NOTE that bbtools reformat.sh can do this
# same task and it's way faster.
#==============================================================================

import sys
import argparse
#import math
import re
from os.path import isfile, expanduser
#from collections import defaultdict
from time import localtime, time
from Basics import messages as ms
from Basics import utils
#import subprocess as sp
import gzip
import hashlib
from os import unlink
from multiprocessing import cpu_count

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

# sam alignment fields
SAM_QNAME = 0
SAM_FLAG = 1
SAM_RNAME = 2
SAM_POS = 3
SAM_MAPQ = 4
SAM_CIGAR = 5
SAM_RNEXT = 6
SAM_PNEXT = 7
SAM_TLEN = 8
SAM_SEQ = 9
SAM_QUAL = 10
SAM_FIRST_ATTR = 11

#==============================================================================
# main
#==============================================================================

def main(args):

	##	
	## check input files
	##
	
	if not isfile(args.fastq):
		ms.error_message("Input file does not exist")
		return 1
	
	
	rres = core(args)

	return 0


def core(args):

	
	sam0 = ["", "4", "*", "0", "0", "*", "*", "0", "0", "", ""]
	linen = 0
	rnum = 0
	tmpname = hashlib.md5(args.fastq).hexdigest()
	samout = "@HD\tVN:1.0\tSO:unsorted\n"
	
	# open input file
	if re.search("\.gz$", args.fastq):
		fin = gzip.open(args.fastq, "r")
	else:
		fin = open(args.fastq, "r")
	
	# open SAM output
	fout = open("{}.sam".format(tmpname), "w")
	fout.write(samout)
	
	t0 = time()
	
	for szl in fin:
		linen += 1
		if linen == 1:
			# read name
			rname = szl.strip()
		elif linen==2:
			# read
			seq = szl.strip()
		elif linen==4:
			qual = szl.strip()
			# get the sam alignment read for writing
			sam = list(sam0)
			sam[SAM_QNAME] = re.sub("^\@", "", rname)
			sam[SAM_SEQ] = seq
			sam[SAM_QUAL] = qual
			# write line out to file
			fout.write("\t".join(sam))
			fout.write("\n")
			# reset line counter
			linen = 0
			rnum += 1
			
		if rnum > 0 and (rnum % 1000000) == 0:
			ms.progress_message("parsed {} reads".format(rnum))
	
	ms.progress_message("parsed {} reads".format(rnum), last=True)
	
	ms.time_diff(t0)
	
	fout.close()
	fin.close()
	
	ms.message("Converting to BAM")
	t0 = time()
	
	t = args.t
	if t==0:
		t = cpu_count()/2
	elif t > cpu_cout():
		t = cpu_count()
	
	cmd = "samtools view -bS -@ {} -o {} {}.sam".format(t, args.bam, tmpname)
	rres =  utils.runcmd(cmd)
	if rres[0] != 0:
		ms.error_message("Failed to create BAM file!")
		return 1
		
	ms.time_diff(t0)
	
	unlink("{}.sam".format(tmpname))
	
	return 0


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Convert FASTQ format alignments to BAM", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('fastq', type=str, help="Input FASTQ")
parser.add_argument('bam', type=str, help="Output BAM file")
parser.add_argument('-t', type=int, default=0, help="BAM compression threads. 0 uses all/2.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

