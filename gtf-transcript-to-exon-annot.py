#!/usr/bin/python
#==============================================================================
# template.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
from Basics import messages as ms
import copy

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
GTF_ATTRIBUTE = 8


#==============================================================================
# main
#==============================================================================

def main(args):

	
	if not isfile(args.gtf):
		ms.error_message("Input file does not exist: {}".format(args.gtf))
		return 1
	
	
	rres = core(args)
	
	

	return rres


def core(args):
	
	dtid = {}
	
	ms.message("Parsing transcript entries in GTF")
	with open(args.gtf, "r") as fin:
		# parse the GTF to find the 'transcript' rows and gather info
		
		for szl in fin:
			if szl[0]=="#":
				continue
			
			grow = gtf_parseline(szl)
			
			if grow['type'] == "transcript":
				tid = grow['attrs']['transcript_id']
				dtid[tid] = grow['attrs']
	
	# finished. now we can pass back through and parse the exon rows out
	szout = ""
	ms.message("Incorporating transcript annotations into exon fields")
	with open(args.gtf, "r") as fin:
		
		for szl in fin:
			if szl[0]=="#":
				szout += szl
				continue
			
			grow = gtf_parseline(szl)
			if grow['type'] == "exon":
				tid = grow['attrs']['transcript_id']
				if tid in dtid:
					# combine annotation information
					rres = merge_annot(grow['attrs'], dtid[tid])
					grow['attrs'] = copy.deepcopy(rres)
				
				szout += gtf_row_tostring(grow)
				szout += "\n"
	
	if args.o is not None:
		with open(args.o, "w") as fout:
			fout.write(szout)
	else:
		sys.stdout.write(szout)
	
	return 0
			
			

def merge_annot(target, source):
	
	for k in source.keys():
		if k not in target:
			target[k] = source[k]
	
	return target

			
def gtf_row_tostring(obj):
	lout = [obj['rname'], obj['db'], obj['type'], str(obj['start']), str(obj['end']), ".", obj['strand'], "."]
	
	# combine the attributes
	attr = obj['attrs']
	lattrs = []
	if "transcript_id" in attr:
		lattrs.append("transcript_id \"{}\"".format(attr['transcript_id']))
	
	if "gene_id" in attr:
		lattrs.append("gene_id \"{}\"".format(attr['gene_id']))
	
	if "gene_name" in attr:
		lattrs.append("gene_name \"{}\"".format(attr['gene_name']))
	
	lnot = set(["transcript_id", "gene_id", "gene_name"])
	
	for aid in attr.keys():
		if aid not in lnot:
			lattrs.append("{} \"{}\"".format(aid, attr[aid]))
		
	lout.append("; ".join(lattrs))
	
	return "\t".join(lout)
	

def gtf_parseline(sz):

	tmp = sz.strip().split("\t")

	grow = {
		'rname':tmp[0],
		'db':tmp[1],
		'type':tmp[2],
		'start':int(tmp[3]),
		'end':int(tmp[4]),
		'strand':tmp[6],
		'attrs':{}}

	# parse attributes
	fsplit = tmp[8].split("\"")
	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		grow['attrs'][key.strip()] = fsplit[i+1].strip()
		i += 2

	return grow

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Take annotation information from the 'transcript' fields in a GTF and incorporate them into the 'exon' fields.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('gtf', type=str, help="GTF annotation to process")
parser.add_argument('-o', type=str, default=None, 
	help="Output file name (results writtin to stdout otherwise")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

