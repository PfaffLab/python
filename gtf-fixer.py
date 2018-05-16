#!/usr/bin/python
#==============================================================================
# gtf-fixer.py
#
# Shawn Driscoll
# 20170614 - happy birthday Mr. President
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys
import argparse
#import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime

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
HBIN = 16000

# GTF fields
GTF_RNAME = 0
GTF_SOURCE = 1
GTF_FEATURE = 2
GTF_START = 3
GTF_END = 4
GTF_SCORE = 5
GTF_STRAND = 6
GTF_FRAME = 7
GTF_ATTR = 8

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	# dict containing the annotation 
	dannot = defaultdict(list)
	# this will be a dict of lists of transcript id and transcript start position so that i can 
	# generate an output that has transcripts ordered by genomic position
	drbins = defaultdict(list)
	lcount = 0
	
	message("Loading {}".format(args.gtf))
	
	# load gtf into a dict keyed by the transcript ids
	with open(args.gtf, "r") as fin:
		for szl in fin:
			lcount += 1
			grow = gtf_parse_line(szl)
			tid = grow[GTF_ATTR]['transcript_id']
			
			# check gene name field
			if "gene_name" in grow[GTF_ATTR]:
				gname0 = grow[GTF_ATTR]['gene_name']
				grow[GTF_ATTR]['gene_name'] = re.sub("[\(\)\{\}\[\]\,\.\s]", "_", gname0)
			
			dannot[tid].append(grow)
			
			if (lcount % 2**14) == 0:
				progress_message("parsed {} lines".format(lcount))
	
	progress_message("parsed {} lines".format(lcount), last=True)
	
	# per transcript id, sort the features by position
	message("sorting exons within transcripts")
	for tid in dannot.keys():
		dannot[tid].sort(key=lambda x: x[GTF_START])
		drbins[dannot[tid][0][GTF_RNAME]].append([tid, dannot[tid][0][GTF_START]])

	# sort the transcripts by position

	lcount = 0
	for rname in sorted(drbins.keys()):
		# dsort the transcripts by position
		drbins[rname].sort(key=lambda x:x[1])
		for ll in drbins[rname]:
			tid = ll[0]
			for grow in dannot[tid]:
				sz = gtf_line_to_string(grow)
				print sz
				lcount += 1
				if (lcount % 2**14) == 0:
					progress_message("wrote {} lines".format(lcount))
		
	progress_message("wrote {} lines".format(lcount), last=True)
			
	return 0



def gtf_parse_line(sz):

	tmp = sz.strip().split("\t")

	for k in [GTF_START, GTF_END]:
		tmp[k] = int(tmp[k])

	# parse attributes
	attrs = {}
	fsplit = tmp[GTF_ATTR].split("\"")
	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	tmp[GTF_ATTR] = attrs
	return tmp

def gtf_line_to_string(grow):
	# collapse the attributes dict down into the string format
	lattr = []
	for k in grow[GTF_ATTR].keys():
		lattr.append("{} \"{}\";".format(k, grow[GTF_ATTR][k]))
	grow[GTF_ATTR] = " ".join(lattr)

	sz = "\t".join(map(str, grow))
	return sz		


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

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Remove unwanted spaces and characters from gene names, reorder exons within transcripts.")
parser.add_argument('gtf', type=str, help="GTF to process")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


