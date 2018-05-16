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
	dtid2gname = {}
	dchrom = defaultdict(list)
	
	fout_name = re.sub("\.gtf$", ".refFlat", args.gtf)
	
	message("Loading {}".format(args.gtf))
	
	# load gtf into a dict keyed by the transcript ids
	with open(args.gtf, "r") as fin:
		for szl in fin:
			if szl[0] == "#":
				continue
			
			grow = gtf_parse_line(szl)
			if grow[GTF_FEATURE] != "exon":
				continue
				
			tid = grow[GTF_ATTR]['transcript_id']
			
			# check gene name field
			if "gene_name" in grow[GTF_ATTR]:
				gname0 = grow[GTF_ATTR]['gene_name']
				grow[GTF_ATTR]['gene_name'] = re.sub("[\(\)\{\}\[\]\,\.\s]", "_", gname0)
			
				dtid2gname[tid] = grow[GTF_ATTR]['gene_name']
			else:
				dtid2gname[tid] = grow[GTF_ATTR]['gene_id']

	# create gene pred format
	message("Creating genePred version of GTF and merging in gene names")
	p1 = sp.Popen("{}/opt/ucsc_tools/gtfToGenePred {} {}.gpred".format(HOME, args.gtf, args.gtf).split())
	p1.wait()
	
	# open that file up and insert the gene names
	message("Writing {}".format(fout_name))
	lpred = []
	idx = 0
	with open("{}.gpred".format(args.gtf), "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			tid = aln[0]
			chrom = aln[1]
			tstart = int(aln[4])
			dchrom[chrom].append([idx, tid, tstart])
			lpred.append(szl)
			idx += 1

	with open(fout_name, "w") as fout:
		for chrom in sorted(dchrom.keys()):
			# sort the transcripts in this chrom
			dchrom[chrom].sort(key=lambda x:x[2])
			for l in dchrom[chrom]:
				tid = l[1]
				gene = dtid2gname[tid]
				fout.write("{}\t{}".format(gene, lpred[l[0]]))
						
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


parser = argparse.ArgumentParser(description="Produce a refFlat annotation from a GTF annotation")
parser.add_argument('gtf', type=str, help="GTF to process")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


