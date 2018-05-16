#!/usr/bin/env python
#==============================================================================
# asi-count-from-junctions.py
#
# Shawn Driscoll
# 20160209
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
#  
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

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
	
	# load the junctions and counts
	sys.stderr.write("Parsing {}\n".format(args.junc_counts))
	djcounts = load_juncs(args.junc_counts, args.b, args.s, args.c)

	# open ref and count stuff
	alid_last = ""
	alid = ""
	lbuffer = []
	fin = open(args.asi_ref, "r")
	
	sys.stderr.write("Parsing asi reference and quantifying\n")
	# print a header out
	print "\t".join(["aloc_id", "path_id", "gene_id", "gene_name", "asi_type", "index_within", "strand", "introns", "transcript_ids", "intron_counts", "total_path_count", "mean_path_count", "aloc_count", "path_ratio"])
	
	for szl in fin:
		aln = szl.strip().split("\t")

		alid = aln[0]
		if alid != alid_last:
			if len(lbuffer) > 0:
				res = process_buffer(lbuffer, djcounts)
				for i in range(len(lbuffer)):
					print "\t".join(map(str, lbuffer[i]))

			lbuffer = []

		alid_last = alid
		lbuffer.append(aln)

	fin.close()

	# process what is left in the buffer
	if len(lbuffer) > 0:
		res = process_buffer(lbuffer, djcounts)
		for i in range(len(lbuffer)):
			print "\t".join(map(str, lbuffer[i]))

	return 0



def load_juncs(f, bed, skip, col):

	djcount = {}
	fin = open(f, "r")
	i = 0
	if skip > 0:
		while i < skip:
			szl = fin.readline()
			i += 1

	for szl in fin:
		aln = szl.strip().split("\t")
		if bed: 
			btmp = map(int, aln[10].split(","))
			aln[1] = int(aln[1])+btmp[0]+1
			aln[2] = int(aln[2])-btmp[1]

		jid = "{}:{}-{}".format(aln[0], aln[1], aln[2])
		
		if bed:
			djcount[jid] = float(aln[4])
		else:
			djcount[jid] = float(aln[col])

	fin.close()
	return djcount

def process_buffer(lbuffer, djc):
	# deal with aloc in buffer
	aloc_count = 0

	for i in range(len(lbuffer)):
		# at each row append count per junction as a comma separated 
		# set, total count, mean count, aloc count, ratio
		ljid = lbuffer[i][7].split("|")
		
		jhits = []
		for jid in ljid:
			# get hit count for this junction
			if jid in djc:
				jhits.append(float(djc[jid]))
			else:
				jhits.append(0)

		jtotal = sum(jhits)
		jmean = jtotal*1.0/len(jhits)

		lbuffer[i] += [",".join(map(str, jhits)), jtotal, jmean, 0, -1]
		aloc_count += jmean

	# now load in the total count for the loc and calc the ratios
	for i in range(len(lbuffer)):
		lbuffer[i][-2] = aloc_count
		if aloc_count > 0:
			lbuffer[i][-1] = lbuffer[i][-3]*1.0/aloc_count
		else:
			lbuffer[i][-1] = "NA"

	return 0

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('asi_ref', type=str, help="ASI reference (from asi-build-index.py)")
parser.add_argument('junc_counts', type=str, 
	help="Junction counts")
parser.add_argument('-b', action="store_const", const=True, default=False, 
	help="Input junctions are in BED format having values for overlap width in column 11")
parser.add_argument('-s', action="store", type=int, default=0, 
	help="Skip this many lines from the top of the junctions file [0]")
parser.add_argument('-c', action="store", type=int, default=3, 
	help="Zero based column index to find the junction counts. Classic bed format would have this in colun 4. [3]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
