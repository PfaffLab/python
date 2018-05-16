#!/usr/bin/env python
#==============================================================================
# asplice-prepare-intron-reference.py
#
# Shawn Driscoll
# 20130809
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Creates a reference of introns to be quantified as intron vs locus from
# junction data (or whatever I guess...as long as it can be compared to a 
# junction in genome coordinates).
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

_NOVEL_PREFIX = "TCONS"

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
	internal_introns = {}
	loc_index = 0
	tid_to_gid = {}
	merge_cuffcompare = False
	additional_annot = False
	aannot = {}

	# check input file
	if not file_exists(args.gtf):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.gtf)
		return 1

	if len(args.t) > 0 or len(args.c) > 0:
		if not file_exists(args.c) or not file_exists(args.t):
			sys.stderr.write("[main] Error: cuffcompare inputs do not exist: %s, %s\n" % (args.c, args.t))
			return 1
		else:
			merge_cuffcompare = True

	# parse the GTF file into a hash by gene_id and within each gene_id
	# the transcript_id

	sys.stderr.write("[main] loading %s\n" % args.gtf)

	fin = open(args.gtf, "r")

	tid = ""
	gid = ""

	for szl in fin: 
		ll = szl.strip().split("\t")

		if ll[2] == "exon":
			attr = parse_gtf_attr(ll[8])
			gid = attr["gene_id"]
			tid = attr["transcript_id"]
			
			# keep a table that can be used to find gene ids for specific transcript ids
			if tid not in tid_to_gid:
				tid_to_gid[tid] = gid
			
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

	if merge_cuffcompare:
		merge_cufflinks_transfrags(gtf, tid_to_gid, args.c, args.t)

	if len(args.a) > 0:
		# load additional annotation file
		additional_annot = True
		aannot = load_additional_annot(args.a)
		if len(aannot.keys()) == 0:
			return 1
	
	
		

	# --
	# loop back through the gtf information by gene id and process
	# each locus.
	# --
	
	sys.stderr.write("[main] processing gene loci\n")
	
	print "aid\tgene_id\tgene_name\tintron\ttype\talternative\tnovel\tannot\tlocus\ttranscript_types\ttranscript_ids\tstrand\tintron_index"
	
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
			itable = {} # table of unique introns
			e_to_i = {} # table to relate exon ids to intron ids
			i_to_e = {} # table connecting intron ids to exon ids

			left_exons = {}
			right_exons = {}
			left_introns = {}
			right_introns = {}
			internal_introns = {}
			
			transcript_introns = {}
			
			transcript_types = {}
			
			strand = "."
			gname = ""
			
			eid_last = ""

			for tid in lg.keys():
				lt = lg[tid]
				
				if len(lt) > 1:
					# this feature has at least one intron
					eid_last = ""
					transcript_introns[tid] = []
					
					transcript_types[tid] = lt[0][1]

					for i in range(len(lt)):
						# -- loop through the transcript and record its exons and introns
						attr = parse_gtf_attr(lt[i][8])
						tid = attr["transcript_id"]
						
						if strand == ".":
							strand = lt[i][6]
						if len(gname) == 0 and "gene_name" in attr:
							gname = attr["gene_name"]
						
						eid = "|".join([lt[i][0], lt[i][3], lt[i][4]])
						
						# append the exon's index to the transcript id
						# etid = tid + ":" + str(i)
						etid = tid
						
						# put exon into the appropriate table						
						if i == 0:
							if eid not in left_exons:
								left_exons[eid] = []
							
							left_exons[eid].append(etid)
						elif i+1 == len(lt):
							if eid not in right_exons:
								right_exons[eid] = []
							
							right_exons[eid].append(etid)

						# put all exons into this main table
						if eid not in etable:
							etable[eid] = []
						
						etable[eid].append(etid)
						
						# make introns
						if i > 0:
							intron = [lt[i][0], int(lt[i-1][4])+1, int(lt[i][3])-1]
							iid = "|".join(map(str, intron))
							
							# -- insert information into the e_to_i and i_to_e tables
							
							if eid_last not in e_to_i:
								e_to_i[eid_last] = set([])
							if eid not in e_to_i:
								e_to_i[eid] = set([])
							
							e_to_i[eid_last].update([iid])
							e_to_i[eid].update([iid])
							
							if iid not in i_to_e:
								i_to_e[iid] = set([])
							
							i_to_e[iid].update([eid, eid_last])
							
							itid = tid + ":" + str(i-1)
							#itid = tid
							
							if i == 1:
								if iid not in left_introns:
									left_introns[iid] = []
								left_introns[iid].append(itid)
								
								if len(lt) == 2:
									# this intron is a start AND an end so add it to both tables
									if iid not in right_introns:
										right_introns[iid] = []
									right_introns[iid].append(itid)
								
							elif i+1 == len(lt):
								if iid not in right_introns:
									right_introns[iid] = []
								right_introns[iid].append(itid)
							else:
								# use these as a control for the end introns. if any of the 
								# introns from the end intron tables are also in this table
								# then they will not be used for alt-end analysis
								if iid not in internal_introns:
									internal_introns[iid] = []
								internal_introns[iid].append(itid)

							# all introns must be added to the itable since end introns
							# can and often will be part of alt cassettes
							if iid not in itable:
								itable[iid] = []
							itable[iid].append(itid)
							
							transcript_introns[tid].append(iid)
						
						eid_last = eid
								
					# finished with transcript
						
			# -- finished gathering and partitioning the exons and introns in this locus
			
			# prune left_introns and right_introns to make sure they are "end" introns that 
			# are not also internal introns
			kill_list = []
			for iid in left_introns:
				if iid in internal_introns:
					kill_list.append(iid)
			
			for i in range(len(kill_list)):
				left_introns.pop(kill_list[i])
			
			kill_list = []
			for iid in right_introns:
				if iid in internal_introns:
					kill_list.append(iid)

			for i in range(len(kill_list)):
				right_introns.pop(kill_list[i])

			ikeys = itable.keys()

			# -- print out rows for this locus
			for iid in sorted(itable.keys()):
				loc_index += 1
				
				# parse out the transcripts this intron is featured in and the corresponding 
				# intron index to each transcript
				tid_list = []
				ttype_list = []
				iindex_list = []
				for ctid in itable[iid]:
					tmp = ctid.split(":")
					tid_list.append(tmp[0])
					iindex_list.append(tmp[1])
					ttype_list.append(transcript_types[tmp[0]])
				
				aid = "AILOC_{:06d}".format(loc_index)
				
				# loose classification - left, right or internal?
				type = ""
				
				if iid in left_introns:
					if strand == "+":
						type = "5p"
					else:
						type = "3p"
				elif iid in right_introns:
					if strand == "+":
						type = "3p"
					else:
						type = "5p"
				else:
					type = "internal"
			
				# determine if this intron shares left or right anchors with any others
				# which identifies this intron as alternative or not and build list of all
				# other introns in the locus
				alternative = False
				other_list = []
				for iid2 in ikeys:
					if iid2 != iid:
						if get_lanchor(iid) == get_lanchor(iid2) or get_ranchor(iid) == get_ranchor(iid2):
							alternative = True
						
						other_list.append(loc_reformat(iid2))
								
				push_out_tab(aid)
				push_out_tab(gid)
				push_out_tab(gname)
				push_out_tab(loc_reformat(iid))
				push_out_tab(type)
				
				if alternative:
					push_out_tab("1")
				else:
					push_out_tab("0")
				
				if merge_cuffcompare:
					# figure out if this is a novel intron
					
					novel_count = 0
					for tid in tid_list:
						if _NOVEL_PREFIX in tid:
							novel_count += 1
					
					if novel_count == len(tid_list):
						push_out_tab("yes")
					else:
						push_out_tab("no") 
											
				else:
					push_out_tab("-")
				
				if additional_annot:
					if iid in aannot:
						push_out_tab(aannot[iid])
					else:
						push_out_tab("-")
				else:
					push_out_tab("-")
				
				push_out_tab(",".join(other_list))
				push_out_tab(",".join(ttype_list))
				push_out_tab(",".join(tid_list))
				push_out_tab(strand)
				push_out_eol(",".join(iindex_list))
				
				#print loc_reformat(iid), get_lanchor(iid), itable[iid], i_to_e[iid]
					
	return 0

def load_additional_annot(fname):
	
	ahash = {}
	iid = ""
	
	if not file_exists(fname):
		sys.stderr.write("[load_additional_annot] Error: failed to open input file\n")
		return {}

	fin = open(fname, "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		iid = "|".join(ll[0:3])
		ahash[iid] = ll[3]
	
	fin.close()
	
	return ahash


def push_out_tab(sz):
	sys.stdout.write(sz + "\t")

def push_out_eol(sz):
	sys.stdout.write(sz + "\n")

def loc_reformat(loc):
	ll = loc.split("|")
	return ll[0] + ":" + ll[1] + "-" + ll[2]

def new_me_alt():
	return {"refs": [], "alts": [], "special": 0, "ref_tid": [], "alt_tid": []}

def new_alt():
	return {"A": [], "B":[], "A_tid": [], "B_tid": [], "special": 0, "novel": 0, "a_novel": False, "b_novel": False}

def get_lanchor(fid):
	feature = fid.split("|")
	return feature[1]

def get_ranchor(fid):
	feature = fid.split("|")
	return feature[2]

def get_edge(fid, idx):
	# 0 is left edge, 1 is right edge
	feature = fid.split("|")
	return int(feature[idx+1])

def feature_length_from_id(fid):
	feature = fid.split("|")
	flen = int(feature[2]) - int(feature[1]) + 1
	return flen

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

def merge_cufflinks_transfrags(rgtf, tid_to_gid, fname_cgtf, fname_ctracking):
	# -- 
	# rgtf is the gid keyed parsed GTF file for the annotation. tid_to_gid is the
	# tid to gene id translation table because we'll be looking at transcript
	# ids from the cuffcompare output and need to get them out of the gene id 
	# keyed gtf information. the final two args are the filenames for the 
	# cuffcompare outputs (typically cuffcmp.combined.gtf and cuffcmp.tracking)
	# -- 
	
	cgtf = {}
	ctracking = {}
	
	# -- parse cufflinks gtf by transcript id
	sys.stderr.write("[merge_cufflinks_transfrags] loading transfrags from %s\n" % fname_cgtf)
	fin = open(fname_cgtf, "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		
		if ll[2] == "exon":
			attr = parse_gtf_attr(ll[8])
			tid = attr["transcript_id"]
			
			if tid not in cgtf:
				cgtf[tid] = []
			
			# add the feature
			if len(cgtf[tid]) > 0:
				# add to front or back of the current list depending on feature start
				i = 0
				while i < len(cgtf[tid]) and int(ll[3]) > int(cgtf[tid][i][3]):
					i += 1
				cgtf[tid].insert(i, list(ll))
			else:
				cgtf[tid].append(list(ll))
	
	fin.close()
	
	# -- 
	# parse tracking. need this simply as a conversion between cufflinks generated
	# transcript ids and those used in the reference gtf. load only 'j' features
	# since those are specifically "novel" configurations that closely match 
	# something annotated.
	# -- 
	sys.stderr.write("[merge_cufflinks_transfrags] processing %s\n" % fname_ctracking)
	
	fin = open(fname_ctracking, "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		ctid = ll[0]
		
		if ll[3] == "j" and ll[2] != "-" and ctid in cgtf:
			rtid = ll[2].split("|")[1]
			rgid = tid_to_gid[rtid]
			n_attr = "gene_id \"%s\"; transcript_id \"%s\";" % (rgid, ctid)
			
			# get this transcript
			c_transcript = cgtf[ctid]
			r_transcript = rgtf[tid_to_gid[rtid]][rtid]
			n_transcript = []
			
			c_first_edge = int(c_transcript[0][4])
			c_last_edge = int(c_transcript[-1][3])
			left_linked = False
			right_linked = False
			left_link_idx = -1
			right_link_idx = -1
			
			if c_first_edge < int(r_transcript[0][4]):
				# first edge of the cufflinks transcript is to the left of the annotated transcript's first 
				# edge

				# find right linkage
				for i in range(len(r_transcript)):
					if c_last_edge == int(r_transcript[i][3]):
						right_linked = True
						right_link_idx = i
						break
				
			elif c_first_edge > int(r_transcript[0][4]):
				# first edge is to the right of the annotated transcript so we can advance until 
				# we find where the cufflinks fragment starts
				
				# figure out where this fragment links in
				
				for i in range(len(r_transcript)):
					if c_first_edge == int(r_transcript[i][4]):
						left_linked = True
						left_link_idx = i
					elif c_last_edge == int(r_transcript[i][3]):
						right_linked = True
						right_link_idx = i
								
			else:
				# aligned with first edge
				left_linked = True
				left_link_idx = 0
				
				# find right linkage
				for i in range(len(r_transcript)):
					if c_last_edge == int(r_transcript[i][3]):
						right_linked = True
						right_link_idx = i
						break
			
			if left_linked or right_linked:
			
				# build the new transcript
				if left_linked:
					# add exons from r_transcript up to and including the one at 
					# left_link_idx
					i = 0
					while i <= left_link_idx:
						n_transcript.append(list(r_transcript[i]))
						n_transcript[-1][8] = n_attr
						i += 1
					
				# now link in the novel transcript from the second exon until one 
				# less than the end
				if len(c_transcript) > 2:
					i = 1
					while i < (len(c_transcript)-1):
						n_transcript.append(list(c_transcript[i]))
						n_transcript[-1][8] = n_attr
						i += 1
				
				# if right_linked then pick up from the annotated transcript at the right_link_idx
				# otherwise include the last exon of the c_transcript and that'll do it
				if right_linked:
					i = right_link_idx
					while i < len(r_transcript):
						n_transcript.append(list(r_transcript[i]))
						n_transcript[-1][8] = n_attr
						i += 1
				else:
					n_transcript.append(list(c_transcript[-1]))
					n_transcript[-1][8] = n_attr
				
				# done with this, now add this transcript into the gene_id bucket so that main can
				# process it.
				rgtf[rgid][ctid] = list(n_transcript)
			
	fin.close()
	
	

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


parser = argparse.ArgumentParser(description="Generates an alternative splicing table from a GTF file to be used for downstream alternative splicing analysis.")
parser.add_argument('gtf', type=str, help="GTF annotation to evaluate relative to the junctions.")
parser.add_argument('-c', type=str, dest="c", default="", help="Cuffcompare 'combined' GTF (forces -t)")
parser.add_argument('-t', type=str, dest="t", default="", help="Cuffcompare 'tracking' file (forces -c)")
parser.add_argument('-a', type=str, dest="a", default="", help="Additional intron annotation in BED format (information is pulled from name column)")
parser.add_argument('--no-novel-vs-novel', dest="no_novel_vs_novel", action="store_const", const=True, default=False, 
				help="Disallow novel vs novel alternative splice locations (default: allowed)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
