#!/usr/bin/env python
#==============================================================================
# exon-useage.py
#
# Shawn Driscoll
# 20130708
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Quantify exon useage relative to gene useage.  
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

_HASH_MOD = 16000

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	elk = {}
	gtf = {}

	# check input file
	if not file_exists(args.alignments):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.alignments)
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

	sys.stderr.write("[main] building unique exon hash\n")
	
	# --
	# at each 'gene_id' cluster, find all unique exons and slam them 
	# into a fast lookup hash
	# --

	for gid in gtf.keys():
		lg = gtf[gid]

		# there has to be at least one multi exon transcript for this to work
		multi_exon_count = 0
		if len(lg) > 0:
			for tid in lg.keys():
				if len(lg[tid]) > 1:
					multi_exon_count += 1


		if multi_exon_count > 0:
			# -- at least one multi-exon transcript, continue...

			# build a table of unique exons
			etable = {}
			for tid in lg.keys():
				lt = lg[tid]

				if len(lt) > 1:
					# this feature as at least one intron

					for i in range(len(lt)):
						eid = "|".join([lt[i][0], lt[i][3], lt[i][4]])

						# append the exon index
						lt[i].append(i)

						if eid not in etable:
							etable[eid] = []

						etable[eid].append(lt[i])
			
			# -- 
			# add these unique exons into the lookup hash
			# --

			# find overlaps of the exons and the introns to find exons that are alternativly 
			# spliced. the search is for exons that are completely contained wihin introns
			ilk = {}
			for iid in itable.keys():
				intron = iid.split("|")
				kid_list = hash_feature(intron)
				
				for kid in kid_list:
					if kid not in ilk:
						ilk[kid] = []

					ilk[kid].append(intron)

			for eid in etable.keys():
				exon = eid.split("|")

				hits = lookup_intersections(exon, ilk)
#				if eid == "13|100131080|100131175":
#					print hits

				if len(hits) > 0:
					# check hits for 100% containment of the exon
					
					for i in range(len(hits)):
						if hits[i][1][1] == 1:

							if eid not in ae_exons:
								ae_exons[eid] = []

							ae_exons[eid].append("|".join(hits[i][0]))

#			if gid == "ENSMUSG00000021645":
#				print etable.keys()
#				print "intron list", len(itable.keys())
#				print itable.keys()
#				print ae_exons
#				return 0

			# if there were alternative exons, continue
			if len(ae_exons.keys()) > 0:
				for eid in ae_exons.keys():

					# find the introns that support this exon
					inc_left = []
					inc_right = []

					exon = eid.split("|")

					# -- find left side junction
					etemp = [exon[0], int(exon[1])-4, int(exon[1])]
					hits = lookup_intersections(etemp, ilk)
					if len(hits) > 0:
						for i in range(len(hits)):
							# check that the intron ends 1 position to the left of this exon
							if int(hits[i][0][2])+1 == int(exon[1]):
								inc_left.append(list(hits[i][0]))

					# -- find right side junction
					etemp = [exon[0], int(exon[2]), int(exon[2])+4]
					hits = lookup_intersections(etemp, ilk)
					if len(hits) > 0:
						for i in range(len(hits)):
							# check that the intron starts 1 position to the right of this exon
							if int(hits[i][0][1])-1 == int(exon[2]):
								inc_right.append(list(hits[i][0]))

#					print "exon", eid
#					print "skips"
#					print ae_exons[eid]
#					print "inc left"
#					print inc_left
#					print "inc right"
#					print inc_right

					# -- find counts for these junctions

					skip_count = 0
					inc_left_count = 0
					inc_right_count = 0
					no_left = len(inc_left) == 0
					no_right = len(inc_right) == 0

					for iid in ae_exons[eid]:
						if iid in juncs:
							for i in range(len(juncs[iid])):
								skip_count += int(juncs[iid][i][3])

					for i in range(len(inc_left)):
						iid = "|".join(inc_left[i])
						if iid in juncs:
							for i in range(len(juncs[iid])):
								inc_left_count += int(juncs[iid][i][3])						

					for i in range(len(inc_right)):
						iid = "|".join(inc_right[i])
						if iid in juncs:
							for i in range(len(juncs[iid])):
								inc_right_count += int(juncs[iid][i][3])						

					# -- build output summary for this exon
					
					# ae id, just the exon location
					lout = []
					lout.append("%s:%s-%s" % (exon[0], exon[1], exon[2]))

					# find gene id and gene name for this event - get from etable with 
					# this exon id
					attr = parse_gtf_attr(etable[eid][0][8])
					lout.append(attr["gene_id"])
					if "gene_name" in attr:
						lout.append(attr["gene_name"])
					else:
						lout.append("-")

					# -- determine type
					strand = etable[eid][0][6]

					if no_left:
						if strand == "+":
							lout.append("5p")
						else:
							lout.append("3p")
					elif no_right:
						if strand == "+":
							lout.append("3p")
						else:
							lout.append("5p")
					else:
						lout.append("internal")

					# -- spacer
					lout.append("-")

					# -- get list of transcript ids this feature is a part of

					tidset = set([])
					for i in range(len(etable[eid])):
						attr = parse_gtf_attr(etable[eid][i][8])
						tidset.update([attr["transcript_id"]])

					lout.append(",".join(list(tidset)))

					# -- make a list of the skipping introns and also assemble a transcript
					# -- list of the transcripts that these introns belong to
					ski = []
					tidset = set([])
					for iid in ae_exons[eid]:
						intron = iid.split("|")
						ski.append("%s:%s-%s" % (intron[0], intron[1], intron[2]))

						for intron in itable[iid]:
							temp_eid = intron[-1]

							for i in range(len(etable[temp_eid])):
								attr = parse_gtf_attr(etable[temp_eid][i][8])
								tidset.update([attr["transcript_id"]])

					lout.append(",".join(ski))
					lout.append(",".join(list(tidset)))

					# -- spacer
					lout.append("-")

					# -- append counts and percent inclusion

					lout.append(inc_left_count)
					lout.append(inc_right_count)
					lout.append(skip_count)

					# include percent inclusion
					inc_final = max(inc_left_count, inc_right_count)
					depth = inc_final + skip_count
					if depth > 0:
						pic = inc_final*1.0/depth
					else:
						pic = 0

					lout.append(pic)

					ae_table[lout[0]] = list(lout)


	# print header
	print "aid\tgene_id\tgene_name\texon_type\texon_index\ttranscript_list\tskipping_intron\tskipping_transcripts\tintron_index\tinc1\tinc2\tskip\tpic"

	# sort and print the ae hits
	for aeid in sorted(ae_table.keys()):
		print "\t".join(map(str, ae_table[aeid]))

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
