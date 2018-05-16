#!/usr/bin/env python
#==============================================================================
# ae-quantify-from-junctions.py
#
# Shawn Driscoll
# 20130708
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Quantify alternative exon useage from junctions and a GTF annotation.  The 
# GTF is used as a way to relate the junctions to exons and to identify
# the alternativly spliced exons from the annotation.  The annotation could
# be generated by something like cufflinks.  The annotation should have a 
# unique gene id for each group of features that should be considered together
# for alternative splicing.  Ensemble GTFs work. 
#==============================================================================

# -- 
# NOTES
# --

# to find the list of exons each junction belongs to within the annotation
# i need to hash the exons and then for each junction make the left anchor
# equal to the start-1 to the start+1 and the right anchor be the end-1 to 
# the end+1. when the overlaps come back find the ones where the exon
# ends at start-1 or the exon starts at end+1.

# i guess i'll find the exons that are skipped first, keep the total skip count
# and then hash those exons to find the junction counts that include them.
# the output will be exon-centric giving the location of the exon, its GTF
# information and the skip/inclusion count information

# many exons will the redundant so their annotations will have to be 
# merged

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
	ae_table = {} # hash of alternative exons
	glk = {}
	jlk = {}
	gtf = {}
	intron_table = {} # table of annotated introns with their hit counts
	juncs = {}

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

	# --
	# parse the junctions into memory so we can snag counts quickly
	# --
	
	sys.stderr.write("[main] parsing junctions\n")

	fin = open(args.juncs, "r")

	for szl in fin:
		ll = szl.strip().split("\t")

		iid = "|".join(ll[0:3])

		if iid in juncs:
			sys.stderr.write("[main] Warning: duplicate junction found: %s\n" % iid)

		juncs[iid] = list(ll)

	fin.close()

	sys.stderr.write("[main] quantifying inclusions for alternative exons\n")


	# --
	# loop back through the parsed information. this time for each gene id
	# group we need to find any alternativly spliced exons. so build a
	# lookup table of the introns, find overlaps with the exons and assemble
	# everything so that we can quantify the intron counts against the junction
	# file supplied.
	# --

	for gid in gtf.keys():
		lg = gtf[gid]

		# for this to work there must be at least one multi-exon isoform
		multi_exon_count = 0
		if len(lg) > 1:
			for tid in lg.keys():
				if len(lg[tid]) > 1:
					multi_exon_count += 1


		if multi_exon_count > 0:
			# continue....

			# build a table of unique exons and introns
			etable = {} # table of unique exons by id containing all GTF information for each exon
			e_to_i = {} # table to relate exon ids to intron ids
			i_to_e = {} # table connecting intron ids to exon ids
			i_hits = {} # table for intron hits
			itable = {} # table of unique introns
			gid_depth = 0 # used as a total of all junction hits for this gene id
			ae_exons = {} # hash of alternative exons relating them to the introns that skip them
			eid_last = ""
			for tid in lg.keys():
				lt = lg[tid]

				if len(lt) > 1:
					# this feature as at least one intron
					eid_last = ""

					for i in range(len(lt)):
						eid = "|".join([lt[i][0], lt[i][3], lt[i][4]])

						# append the exon index
						lt[i].append(i)

						if eid not in etable:
							etable[eid] = []

						etable[eid].append(lt[i])

						if i > 0:
							# make the intron feature - it's position and the ids of both exons
							# it connects
							# intron = [lt[i][0], int(lt[i-1][4])+1, int(lt[i][3])-1, eid_last, eid, 0]
							iid = "|".join(map(str, [lt[i][0], int(lt[i-1][4])+1, int(lt[i][3])-1]))
							
							if iid not in i_hits:
								# insert iid
								i_hits[iid] = 0
								# grab count from juncs if its there
								if iid in juncs:
									i_hits[iid] = int(juncs[iid][3])
									gid_depth += int(juncs[iid][3])

							# create link to exons and this junction in e_to_i table for quick access later
							
							if eid_last not in e_to_i:
								e_to_i[eid_last] = set([])
							
							if eid not in e_to_i:
								e_to_i[eid] = set([])
							
							if iid not in i_to_e:
								i_to_e[iid] = set([])
							
							e_to_i[eid_last].update([iid])
							e_to_i[eid].update([iid])
							i_to_e[iid].update([eid_last, eid])

						eid_last = eid
			
			
#			if gid == "ENSMUSG00000020140":
#				for iid in i_hits:
#					print iid, i_hits[iid]
				
#				print gid_depth
#				return 1
			
			# --
			# make a fast lookup hash of the introns ... probably overkill but whatever
			# --
			ilk = {}
			for iid in i_hits.keys():
				intron = iid.split("|")
				kid_list = hash_feature(intron)
				
				for kid in kid_list:
					if kid not in ilk:
						ilk[kid] = []

					ilk[kid].append(intron)
			
			# --
			# determine which exons are contained by junctions - these are the ones that we'll 
			# quantify for testing
			# --
			for eid in etable.keys():
				exon = eid.split("|")

				hits = lookup_intersections(exon, ilk)

				if len(hits) > 0:
					# check hits for 100% containment of the exon
					
					for i in range(len(hits)):
						if hits[i][1][1] == 1:

							if eid not in ae_exons:
								ae_exons[eid] = []
							
							# append intron id for this exon id
							ae_exons[eid].append("|".join(hits[i][0]))

#			if gid == "ENSMUSG00000021645":
#				print etable.keys()
#				print "intron list", len(itable.keys())
#				print itable.keys()
#				print ae_exons
#				return 0

#			if gid == "ENSMUSG00000006732":
#				for eid in ae_exons:
#					print "#----------------- " + eid
#					print etable[eid]
			
#				return 1

			# --
			# check the discovered skipped exons and their intron lists. keep only
			# exons with introns that link back into the same transcript. in other words
			# at least one of the anchors of the skipping intron must also exist
			# in the transcript of the skipped exon
			# -- 
			
			ae_exons_temp = {}
			if len(ae_exons.keys()) > 0:
				for eid in ae_exons.keys():
					
					# get the introns that include this exon
					iset = list(e_to_i[eid])
					
					# make a set out of the starts and ends of the introns
					anchor_set = set([])
					for iid in iset:
						intron = iid.split("|")
						anchor_set.update(intron[1:3])
					
					# loop through skipping introns and see if there is an overlap of 
					# anchors
					for iid in ae_exons[eid]:
						intron = iid.split("|")
						
						rres = anchor_set.intersection(set(intron[1:3]))
						if len(rres) > 0:
							if eid not in ae_exons_temp:
								ae_exons_temp[eid] = []
							
							ae_exons_temp[eid].append(iid)

			
#			if gid == "ENSMUSG00000021645":	
#				for eid in ae_exons:
#					print "#----------- " + eid
#					print ae_exons[eid]
#				print "\n"
#				for eid in ae_exons_temp:
#					print "#----------- " + eid
#					print ae_exons_temp[eid]
							
#				return 1

			# replace alternative exon list with the revised one
			ae_exons = ae_exons_temp
			
			# if there were alternative exons, continue
			if len(ae_exons.keys()) > 0:
				for eid in ae_exons.keys():
					
					exon = eid.split("|")
					
					# for each eid accumulate the inclusion count from its supporting introns from
					# the e_to_i table

#					if gid == "ENSMUSG00000020140" and eid == "10|115478387|115478602":
#						print e_to_i[eid]
#						return 1				

					inclusion_count = 0
					inclusion_num = 0
					skip_count = 0
					
					for iid in list(e_to_i[eid]):
						if i_hits[iid] > 0:
							inclusion_num += 1
						inclusion_count += i_hits[iid]
					
					# inclusion can be the mean of the non-zero intron hits
					if inclusion_num > 0:
						inclusion_count = int(round(inclusion_count*1.0/inclusion_num))
					
					# gather the skip count
					for iid in ae_exons[eid]:
						skip_count += i_hits[iid]

					# -- 
					# get list of transcript ids this feature is a part of and for each
					# transcript figure out if the exon is a 3p, 5p or internal
					# --

					tidset = set([])
					typeset = set([])
					for i in range(len(etable[eid])):
						attr = parse_gtf_attr(etable[eid][i][8])
						tid = attr["transcript_id"]
						tidset.update([tid])
						
						# find number of exons in this transcript
						num_exons = len(gtf[gid][tid])
						strand = etable[eid][i][6]
						eidx = int(etable[eid][i][-1])+1

						if strand == "+" and eidx == 1:
							typeset.update(["5p"])
						elif strand == "+" and eidx == num_exons:
							typeset.update(["3p"])
						elif strand == "-" and eidx == 1:
							typeset.update(["3p"])
						elif strand == "-" and eidx == num_exons:
							typeset.update(["5p"])
						else:
							typeset.update(["internal"])
											
					typeset = list(typeset)
					tidset = list(tidset)

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

					# -- append type
					lout.append(",".join(typeset))

					# -- spacer
					lout.append("-")

					# -- transcript ids for this exon
					lout.append(",".join(tidset))

					# -- make a list of the skipping introns and also assemble a transcript
					# -- list of the transcripts that these introns belong to
					ski = []
#					tidset = set([])
					for iid in ae_exons[eid]:
						intron = iid.split("|")
						ski.append("%s:%s-%s" % (intron[0], intron[1], intron[2]))

#						for intron in itable[iid]:
#							temp_eid = intron[-1]

#							for i in range(len(etable[temp_eid])):
#								attr = parse_gtf_attr(etable[temp_eid][i][8])
#								tidset.update([attr["transcript_id"]])

					lout.append(",".join(ski))
#					lout.append(",".join(list(tidset)))
#					lout.append("-")
					lout.append("-")

					# -- spacer
					lout.append("-")

					# -- append counts and percent inclusion

					lout.append(0)
					lout.append(inclusion_count)
					lout.append(skip_count)

					# include percent inclusion
#					inc_final = max(inc_left_count, inc_right_count)
#					depth = inc_final + skip_count
#					if depth > 0:
#						pic = inc_final*1.0/depth
#					else:
#						pic = 0	

#					if gid_depth > 0:
#						pic = inclusion_count*1.0/gid_depth
#					else:
#						pic = 0

					depth = inclusion_count + skip_count
					if depth > 0:
						pic = inclusion_count*1.0/depth
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
