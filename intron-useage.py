#!/usr/bin/env python
#==============================================================================
# intron-useage.py
#
# Shawn Driscoll
# 20130708
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Quantify intron useage (intron hits vs total hits to introns for a single
# gene). 
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

# -- 
# globals
# --

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	gtf = {}
	juncs = {}
	
	gid_hits = 0
	

	# check input file
	if not file_exists(args.juncs):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.juncs)
		return 1
	if not file_exists(args.gtf):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.gtf)
		return 1


	# parse the GTF file into a hash by gene_id and within each gene_id
	# the transcript_id

	sys.stderr.write("[main] parsing GTF\n")

	fin = open(args.gtf, "r")

	tid = ""
	gid = ""

	for szl in fin: 
		ll = szl.strip().split("\t")

		if ll[2] == "exon":
			attr = parse_gtf_attr(ll[8])
			gid = attr["gene_id"]
			tid = attr["transcript_id"]

			# add gene id if it's not in the table
			if gid not in gtf:
				gtf[gid] = {}

			# add the transcript id if its not in the gene id table
			if tid not in gtf[gid]:
				gtf[gid][tid] = []

			# add the feature
			if len(gtf[gid][tid]) > 0:
				# add to front or back of the current list
				i = 0
				while i < len(gtf[gid][tid]) and int(ll[3]) > int(gtf[gid][tid][i][3]):
					i += 1
				gtf[gid][tid].insert(i, list(ll))
			else:
				gtf[gid][tid].append(list(ll))

	fin.close()

	sys.stderr.write("[main] parsing junctions\n")

	fin = open(args.juncs, "r")

	for szl in fin:
		ll = szl.strip().split("\t")

		iid = "|".join(ll[0:3])

		if iid not in juncs:
			juncs[iid] = []

		juncs[iid].append(list(ll))

	fin.close()

	sys.stderr.write("[main] quantifying hits at all unique introns\n")


	# --
	# loop back through the parsed information. this time for each gene id
	# group we need to find any alternativly spliced exons. so build a
	# lookup table of the introns, find overlaps with the exons and assemble
	# everything so that we can quantify the intron counts against the junction
	# file supplied.
	# --

	for gid in gtf.keys():
		lg = gtf[gid]
		gid_hits = 0
		gene_name = set([])
		
		# feature needs to have at least one isoform with more than one
		# intron.
		
		multi_exon_count = 0
		
		if len(lg) > 0:
			for tid in lg.keys():
				if len(lg[tid]) > 1:
					multi_exon_count += 1


		if multi_exon_count > 1:
			# -- this gene has at least two isoforms with more than one exon

			# build a table of unique introns
			itable = {}
			for tid in lg.keys():
				lt = lg[tid]
				attr = parse_gtf_attr(lt[0][8])
				if "gene_name" in attr:
					gene_name.update([attr["gene_name"]])
					
				
				if len(lt) > 1:
					# this feature has at least one intron, loop through it
					for i in range(1, len(lt)):

						# make the intron feature - it's position and the ids of both exons
						# it connects
						intron = [lt[i][0], int(lt[i-1][4])+1, int(lt[i][3])-1]
						iid = "|".join(map(str, intron[0:3]))

						if iid not in itable:
							itable[iid] = 0
			
			# check length of the itable, if it is only 1 then skip this gene
			if len(itable) > 1:
				# look up the hits for these introns
				for iid in itable:
					if iid in juncs:
						for i in range(len(juncs[iid])):
							itable[iid] += int(juncs[iid][i][3])
							gid_hits += int(juncs[iid][i][3])
				
				# now all hits are tally'd up, generate simple output table
				# intron, gene_id, gene_name, intron_hits, gene_hits, p_inc
				for iid in itable:
					lout = []
					
					temp = iid.split("|")
					iid_loc = "%s:%s-%s" % (temp[0], temp[1], temp[2])
					
					lout.append(iid_loc)
					lout.append(gid)
					
					if len(gene_name) > 0:
						lout.append(",".join(list(gene_name)))
					else:
						lout.append("-")
					
					lout.append(itable[iid])
					lout.append(gid_hits)
					
					if gid_hits > 0:
						lout.append(itable[iid]*1.0/gid_hits)
					else:
						lout.append(0)
					
					print "\t".join(map(str, lout))


	return 0

def parse_gtf_attr(field):

	attrs = {}

	# split into fields
	l1 = field.split(";")

	# split each field into key and value
	l2 = []
	for i in range(len(l1)):
		l2.append(l1[i].split("\""))


	for i in range(len(l2)):
		if len(l2[i]) > 1:
			attrs[l2[i][0].strip()] = l2[i][1].strip()

	return(attrs)		


def hash_feature(lf):
	# --
	# considers the possibility that a feature would hash into several bins,
	# which totally happens with introns because sometimes they are long
	# -- 

	mstart = int(lf[1])/_HASH_MOD
	mend = int(lf[2])/_HASH_MOD

	kid_list = []

	kid_list.append(lf[0]+str(mstart))
	
	i = mstart+1
	while i <= mend:
		kid_list.append(lf[0]+str(i))
		i += 1

	return kid_list


def hash_pos(rname, pos):

	bucket = int(pos)/_HASH_MOD
	kid = rname + str(bucket)
	return kid

def lookup_intersections(feature, lk):
	# --
	# find intersections of feature (expected to be a list with ref, start, and end)
	# with features in the lookup hash, lk
	# --

	res = []

	feature_set = set([])

	# hash the start and end positions of the feature
#	kid1 = hash_pos(feature[0], feature[1])
#	kid2 = hash_pos(feature[0], feature[2])

	kid_list = hash_feature(feature)

	for kid in kid_list:
		# find overlaps, if any
		if kid in lk:
			lk_features = lk[kid]
			for i in range(len(lk_features)):
				if feature_overlap(feature, lk_features[i][0:3]):
					# make sure we don't double add a feature that's in 
					# more than one bin of the lookup hash
					fid = "".join(lk_features[i][0:3])
					if fid not in feature_set:
						feature_set.update([fid])

						# -- features overlap, find overlap length
						ovlen = overlap_length(feature, lk_features[i][0:3])
						lhit = [list(lk_features[i]), ovlen]
						res.append(lhit)
	

	return res


def feature_overlap(f1, f2):

	if f1[0] == f2[0]:
		# -- same reference name
		if int(f1[1]) < int(f2[2]) and int(f1[2]) > int(f2[1]):
			return True

	return False

def overlap_length(f1, f2):
	# --
	# assumes the features overlap and only seeks to obtain the 
	# length of the overlapped region. returns a list including
	# the length of the overlap and the ratio of the overlap 
	# relative to the lengths of f1 and f2
	# --

	# -- variables
	olen = 0

	# -- we have two features f1 = a, b; f2 = c, d
	# -- overlap condition is a <= d && b >= c

	a = int(f1[1])
	b = int(f1[2])
	c = int(f2[1])
	d = int(f2[2])
	ab = b-a+1
	cd = d-c+1

	if a == c and b == d:
		# -- features are identical
		olen = ab
	else:
		# -- features are not identical

		if a <= c and b >= d:
			# -- cd is contained by ab
			olen = cd
		elif c <= a and d >= b:
			# -- ab is contained by cd
			olen = ab
		else:
			# overlap is partial
			d1 = d-a+1
			d2 = b-c+1
			olen = min(d1, d2)

	return [olen, olen*1.0/ab, olen*1.0/cd]

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


parser = argparse.ArgumentParser(description="Evaluates alternative exon useage directly from junctions.")
parser.add_argument('juncs', type=str, help="Junctions file (from bam-to-junctions)")
parser.add_argument('gtf', type=str, help="GTF annotation to evaluate relative to the junctions.")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
