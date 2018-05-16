#!/usr/bin/env python
#==============================================================================
# asplice-count-hits.py
#
# Shawn Driscoll
# 20130809
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# used to count hits to introns and loci from a .junc file (from bam-to-junctions)
# against a aintron file from asplice-prepare-intron-reference.py 
#==============================================================================

import sys, argparse

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	i_hits= {}
	
	args.m = True

	# check input file
	if not file_exists(args.juncs):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.juncs)
		return 1

	if not file_exists(args.ref):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.ref)
		return 1

	# -- load the junctions
	sys.stderr.write("[main] loading junctions from %s\n" % args.juncs)
	fin = open(args.juncs, "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		
		iid = intron_to_string(ll[0:3])
		if iid in i_hits:
			sys.stderr.write("[main] Error: repeated junction in junction's file: %s\n" % iid)
			return 1
		
		# store count 
		i_hits[iid] = int(ll[3])
		
	fin.close()
	
	
	# -- read through the asplice file and print out report of hits as we go
	sys.stderr.write("[main] quantifying\n")
	
	print "aid\tgene_id\tgene_name\tintron\ttype\talternative\tnovel\tannot\ttranscript_types\ttranscript_ids\tstrand\tintron_index\tintron_hits\tlocus_hits"
	
	fin = open(args.ref, "r")
	szl = fin.readline()
	for szl in fin:
		ll = szl.strip().split("\t")
		
		if (args.a and ll[5] == "1") or not args.a:

			ref_iid = ll[3]
			locus_iid = ll[8].split(",")
					
			ref_hits = 0
			locus_hits = 0
			
			if ref_iid in i_hits:
				ref_hits = i_hits[ref_iid]
			
			for iid in locus_iid:
				if iid in i_hits:
					locus_hits += i_hits[iid]
	
			# drop the locus list
			ll.pop(8)
			
			ll.append(ref_hits)
			ll.append(locus_hits)
			
			print "\t".join(map(str, ll))

	fin.close()

	return 0

def intron_to_string(li):
	return "%s:%s-%s" % (li[0], li[1], li[2])

def geo_mean(x):
	n = len(x)
	root = 1./n
	
	p = 1
	for i in range(n):
		p = p*(x[i]+1)
	
	m = p**root - 1
	return m
	

def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('ref', type=str, help="Alternative splice reference from asplice-prepare-reference")
parser.add_argument('juncs', type=str, help="Junctions file from bam-to-junctions (.junc)")
parser.add_argument('-a', dest="a", action="store_const", const=True, default=False, help="Quantify alternative introns only (default: all)")
args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
