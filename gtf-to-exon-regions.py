#!/usr/bin/python
#==============================================================================
# gtf-to-exon-regions.py
#
# Shawn Driscoll
# 20170601
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse GTF and create a BED format file with exon regions. 
# This will join overlapping regions into single rows.  
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

GTF_RNAME = 0
GTF_SOURCE = 1
GTF_FEATURE = 2
GTF_START = 3
GTF_END = 4
GTF_SCORE = 5
GTF_STRAND = 6
GTF_FRAME = 7
GTF_ATTRIBUTE = 8

# region fields
REGION_RNAME = 0
REGION_START = 1
REGION_END = 2
REGION_TAG = 3

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	dref = None
	foutname = "{}.bed".format(args.o)
	
	#
	# load features
	message("Loading {}".format(args.gtf))
	dref = parse_gtf_to_chrom_lists(args.gtf, args.f)
	# sort chroms
	message("Sorting features per chrom by start position")
	sort_chroms(dref)

	message("Collapsing regions and writing to {}".format(foutname))
	with open(foutname, "w") as fout:
		for k in sorted(dref.keys()):
			message("processing {}".format(k))
			lhat = collapse_chrom_regions(dref[k])
			for r in lhat:
				fout.write(region_tostring(r, use_pid=args.uid)+"\n")	
	

	return 0


#==============================================================================
# defs
#==============================================================================

def print_region_list(l, use_pid=False):
	n = len(l)
	
	for r in l:
		print region_tostring(r, use_pid=use_pid)
	
	return 0

def region_init(rname, start, end, tag=None):
	# region is just a list
	r = [rname, int(start), int(end), None]
	# set optional fields
	if tag is not None:
		r[REGION_TAG] = tag
	
	return r

def region_pid(r):
	return "{}:{}-{}".format(r[REGION_RNAME], r[REGION_START], r[REGION_END])

def region_tolist(r, use_pid=False):
	lout = [r[REGION_RNAME], r[REGION_START], r[REGION_END]]
	
	if use_pid:
		lout.append(region_pid(r))
		
	elif r[REGION_TAG] is not None:
		lout.append(",".join(list(r[REGION_TAG])))
	
	return lout

def region_tostring(r, use_pid=False):
	l = region_tolist(r, use_pid=use_pid)
	return "\t".join(map(str, l))

#
# compare two regions. 
# return value:
# 0 for no overlap
# 1 for overlap
# 2 for identical
def compare_regions(r1, r2):

	rres = 0

	# check ref names. if not equal then we're done
	if r1[REGION_RNAME] != r2[REGION_RNAME]:
		return 0

	# ref names must be equal
	if r1[REGION_START]==r2[REGION_START] and r1[REGION_END]==r2[REGION_END]:
		# starts and ends are identical
		return 2

	# now check for overlap
	if r1[REGION_END] >= r2[REGION_START] and r2[REGION_END] >= r1[REGION_START]:
		# overlap!
		return 1
	
	return rres
	
#
# parse a gtf and build lists per chromosome of all features of type, ftype
def parse_gtf_to_chrom_lists(fname, ftype):
	
	# open fname and load it up
	dchroms = defaultdict(list)
	
	with open(fname, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			if aln[GTF_FEATURE] == ftype:
				# keep it
				rres = re.search("transcript_id \"([^\"]+)\"", szl)
				if rres:
					tid = rres.group(1)
				else:
					tid = "NOID"
					
				r = region_init(aln[GTF_RNAME], aln[GTF_START], aln[GTF_END], tag=set([tid]))
				dchroms[aln[GTF_RNAME]].append(r)
	
	return dchroms

#
# for sorting
def region_sort_key(item):
	return item[REGION_START]

# 
# provided 'd' is the output of parse_gtf_to_chrom_lists			
def sort_chroms(d):
	
	k = None
	
	for k in d.keys():
		l0 = list(d[k])
		# sort it
		d[k] = sorted(l0, key=region_sort_key)
	
	return 0

def collapse_chrom_regions(l):
	lhat = []
	n = len(l)
	
	r = None
	
	for i in range(n):
		# get current region
		r0 = l[i]
		if r is None:
			# set current to the current running region
			r = list(l[i])
		else:
			# we have a running region. see if the one from the chrom list
			# overlaps the running region
			rres = compare_regions(r0, r)
			if rres > 0:
				# yes they overlap. update 'r'
				r[REGION_START] = min([r[REGION_START], r0[REGION_START]])
				r[REGION_END] = min([r[REGION_END], r0[REGION_END]])
				r[REGION_TAG].update(r0[REGION_TAG])
			else:
				# they do not overlap. the running region is finished and we 
				# can append it to the new list
				lhat.append(r)
				r = None
	
	return lhat

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Parse GTF and create a BED format file with exon regions. This will join overlapping regions into single rows.")
parser.add_argument('gtf', type=str, help="GTF file to process")
parser.add_argument('-f', type=str, default="exon", 
	help="Feature type to use for creating regions [exon]")
parser.add_argument('-o', type=str, default="collapsed_gtf", 
	help="Output stub for collapsed bed file")
parser.add_argument('--uid', action="store_const", const=True, default=False, 
	help="Generate per-region unique IDs. Otherwise IDs will be concatenated transcript id values")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
#	except KeyboardInterrupt:
#		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print e
		sys.exit(1)

