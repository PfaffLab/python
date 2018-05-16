#!/usr/bin/env python
#==============================================================================
# juncdb-to-tchains.py
#
# Shawn Driscoll
# 20160330
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses juncdb juncdb quantification output and reformats into a transcript
# keyed table showing the intron chain counts.
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
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
	junc = []
	lres = []
	tlist = []
	lcounts = []
	ljuncs = []
	complete = ""
	thash = {}
	ghash = {}
	jid = ""
	gid = ""
	thash_unique_count = {}

	# parse juncdb file
	sys.stderr.write("parsing input file {}...\n".format(args.juncdb))
	thash, ghash = parse_juncdb(args.juncdb)
	sys.stderr.write("done.\n")

	sys.stderr.write("counting unique introns for transcripts...\n")
	# process the gid hash to make transcript keyed dict indicating their number of unique junctions
	for gid in ghash.keys():
		for jid in ghash[gid].keys():
			if len(ghash[gid][jid])==1:
				tid = ghash[gid][jid][0]
				if tid not in thash_unique_count:
					thash_unique_count[tid] = 0

				thash_unique_count[tid] += 1
	sys.stderr.write("done.\n")

	# print a header
	print "\t".join(["transcript_id", "gene_id", "gene_name", "strand", "total_juncs", "total_hits", "junctions_hit", "complete_ratio",
		"unique_juncs", "unique_reads", "unique_juncs_hit", "ratio_unique_juncs_hit",
		"jmapp_wmean", "tmapp_a", "tmapp_b", "hit_range", "complete", "count_chain", "weight_chain", "junc_chain"])

	sys.stderr.write("generating final output...\n")

	# build output
	for tid in sorted(thash.keys()):
		tlist = sort_transcript(thash[tid])
		
		# loop through and built list of the junctions and a list of their counts
		ljuncs = []
		ljweights = []
		lcounts = []
		num_hit = 0
		total_hits = 0
		total_hits_unique = 0
		total_unique = 0
		total_unique_hit = 0
		total_weight = 0
		total_wcount = 0
		weighted_counts = []


		min_hits = 1e16
		max_hits = 0
		complete = "no"
		tid_unique = {}

		for junc in tlist:
			# make junction id
			jid = "{}:{}-{}".format(junc[0], junc[1], junc[2])
			ljuncs.append(jid)

			# gene gene id
			gid = junc[6].split(",")[0]

			# get weight for this junction

			ljweights.append(1.0/len(ghash[gid][jid]))

			k = float(junc[3])
			lcounts.append(k)
			weighted_counts.append(k*ljweights[-1])
			total_wcount += weighted_counts[-1]
			total_weight += ljweights[-1]

			if ljweights[-1]==1:
				total_unique += 1
				if k > 0:
					total_unique_hit += 1
					total_hits_unique += k

			if k > 0:
				num_hit += 1
				total_hits += k

			if k < min_hits:
				min_hits = k
			if k > max_hits:
				max_hits = k

		
		sz_counts = ",".join(map(str, lcounts))
		sz_juncs = ",".join(ljuncs)
		sz_weights = ",".join(vnumeric_to_string(ljweights))
		ratio_hit = num_hit*1.0/len(lcounts)
		sz_range = ",".join(vnumeric_to_string([min_hits, max_hits]))
		total_juncs = len(ljuncs)
		mweighted_mean_depth = sum(weighted_counts)/total_weight

		tmap_a = total_weight/len(ljuncs)
		tmap_b = 0
		if tid in thash_unique_count:
			tmap_b = thash_unique_count[tid]*1.0/len(ljuncs)

		# is it complete?
		if ratio_hit==1:
			complete = "yes"

		unique_hit_ratio = 0
		if total_unique==0:
			total_unique = "NA"
			total_hits_unique = "NA"
			total_unique_hit = "NA"
			unique_hit_ratio = "NA"
		else:
			unique_hit_ratio = fformat(total_unique_hit*1.0/total_unique)

		# build a list for printing this transcript out
		# tid, gid, gname, strand
		lres = [tid, tlist[0][6], tlist[0][5], tlist[0][4], 
			total_juncs,
			total_hits, 
			num_hit,
			ratio_hit,
			total_unique,
			total_hits_unique,
			total_unique_hit,
			unique_hit_ratio,
			fformat(mweighted_mean_depth),
			fformat(tmap_a),
			fformat(tmap_b),
			sz_range,
			complete,
			sz_counts,
			sz_weights,
			sz_juncs]

		print "\t".join(map(str, lres))

	sys.stderr.write("done.\n")

	return 0


#==============================================================================
# functions
#==============================================================================

def parse_juncdb(f):

	# variables
	aln = []
	szl = ""
	tids = []
	tid = ""
	gid = ""
	thash = {}
	ghash = {}

	# open file
	fin = open(f, "r")
	# skip header
	aln = fin.readline().strip().split("\t")
	tid_idx = 0
	for k in aln:
		if k == "transcript_id":
			break
		tid_idx += 1


	# read it in. we are going to make lists of lists for each transcript id
	# found in the transcript id column

	# we also need a dict by gene id that has a dict of junctions and their associated transcript ids.
	# this will be used to figure out a "mappability" kind of ratio and to track junctions that are 
	# specific to transcripts

	for szl in fin:
		# strip and split line by tabs
		aln = szl.strip().split("\t")

		# split the transcript id field
		tids = aln[tid_idx].split(",")

		# get the gene id
		lgid = aln[tid_idx-1].split(",")

		# add slot for the gene id if necessary
		for gid in lgid:
			if gid not in ghash:
				ghash[gid] = {}

		# make a junction key for this row
		jkey = "{}:{}-{}".format(aln[0], aln[1], aln[2])
		# put the transcript ids from this row into this slot
		for gid in lgid:
			ghash[gid][jkey] = list(tids)

		# loop through transcript ids in this row and create  slots for them. populatte
		# with the junction info from this row
		for tid in tids:
			if tid not in thash:
				thash[tid] = []

			# copy position, count and strand, gene name and gene id
			ltmp = list(aln[0:4] + [aln[5]] + aln[(tid_idx-2):tid_idx])
			thash[tid].append(ltmp)

	fin.close()

	return thash, ghash

# --
# this function sorts a transcript list for a single transcript paying
# attention to strand. if strand is negative then the sort is reversed
def sort_transcript(ltid):

	lout = []
	o = []
	strandPos = True
	v = []

	# first get the start positions from each junction
	for ll in ltid:
		v.append(int(ll[1]))

	# get sort order by coordinate
	o = order(v)

	# check strand
	if ltid[0][4] == "-":
		# reverse the order
		o.reverse()

	# populate a new sorted list
	for i in o:
		lout.append(list(ltid[i]))

	# return it
	return(lout)


# like the R order function which returns a sorting of the indices
def order(v):
	return list(np.argsort(v))


def vnumeric_to_string(v, prec=2):
	lout = []
	fstring = "{:0." + str(prec) + "f}"
	for i in range(len(v)):
		lout.append(fstring.format(v[i]))
	return(lout)

def fformat(n, prec=4):
	fstring = "{:0." + str(prec) + "f}"
	return(fstring.format(n))

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('juncdb', type=str, 
	help="juncdb quantification (parse-run juncdb)")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
