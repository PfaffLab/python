#!/usr/bin/env python
#==============================================================================
# asplice-prepare-reference.py
#
# Shawn Driscoll
# 20130716
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Processes a gene annotation format (GTF) and assembles alternative splice
# arrangements from the transcript models within each gene locus. so far
# it does alt starts and ends (3' and 5' exons), alt cassettes and 
# mutually exclusive exons (in two categories). The output table includes
# location ids and lists of introns that may be used to capture hits
# from junctions or by some other quantification method.
# 20130718
# Added in functionality to import cuffcompare output (compared to the
# same GTF being analyzed by this script) into the analysis.  Fragments
# are integrated into the annotated transcripts they are similar to
# based on cuffcompare's "tracking" file.
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

	# --
	# loop back through the gtf information by gene id and process
	# each locus.
	# --
	
	sys.stderr.write("[main] processing gene loci\n")
	
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
			
			strand = "."
			gname = ""
			
			eid_last = ""

			for tid in lg.keys():
				lt = lg[tid]
				
				if len(lt) > 1:
					# this feature has at least one intron
					eid_last = ""
					transcript_introns[tid] = []

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
							
							#itid = tid + ":" + str(i-1)
							itid = tid
							
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

			# -- assemble alternative splicing events
			
			left_alts = find_alt_ends(left_introns, itable, left_exons, e_to_i, i_to_e, "left")
			right_alts = find_alt_ends(right_introns, itable, right_exons, e_to_i, i_to_e, "right")
			alt_cas = find_alt_cassettes(itable, transcript_introns)
			me_pairs = find_me_pairs(itable, etable, i_to_e, transcript_introns)
			retained_introns = find_retained_introns(itable, etable, i_to_e)
			
			# make sure none of the same introns are used as references in the left and
			# right side alts as in the alt_cas
			cassette_refs = set([])
			for i in range(len(alt_cas)):
				cassette_refs.update(alt_cas[i]["B"])
			
			i = len(left_alts)
			while i:
				i -= 1
				if left_alts[i]["A"][0] in cassette_refs:
					left_alts.pop(i)
			
			i = len(right_alts)
			while i:
				i -= 1
				if right_alts[i]["A"][0] in cassette_refs:
					right_alts.pop(i)
			
#			print "#------------- left_alts"
#			print left_alts
#			print "#------------- right_alts"
#			print right_alts
#			print "#------------- alt_cassettes"
#			print alt_cas
#			print "#------------- me_pairs"
#			print me_pairs
#			return 1

			# print out an annotation for these features
			for i in range(len(left_alts)):
				
				print_me = True
				if args.no_novel_vs_novel:
					if left_alts[i]["a_novel"] == True and left_alts[i]["b_novel"] == True:
						print_me = False
				
				if print_me:
				
					aloc = "ALOC_{:08d}".format(loc_index)
					loc_index += 1
					lout = [aloc, gid, gname, strand]
					
					if left_alts[i]["special"] == 1:
						if strand == "-":
							lout.append("alt_3p")
						else:
							lout.append("alt_5p")
					else:
						if strand == "-":
							lout.append("alt_end")
						else:
							lout.append("alt_start")
					
					if left_alts[i]["novel"] == 1:
						lout[-1] += "*"
					
					lout.append(",".join(left_alts[i]["A_tid"]))
					lout.append(",".join(left_alts[i]["B_tid"]))
					
					lout.append(",".join(map(loc_reformat, left_alts[i]["A"])))
					lout.append(",".join(map(loc_reformat, left_alts[i]["B"])))
					
					print "\t".join(lout)

			for i in range(len(right_alts)):
				
				print_me = True
				if args.no_novel_vs_novel:
					if right_alts[i]["a_novel"] == True and right_alts[i]["b_novel"] == True:
						print_me = False
				
				if print_me:
					aloc = "ALOC_{:08d}".format(loc_index)
					loc_index += 1
					lout = [aloc, gid, gname, strand]
					
					if right_alts[i]["special"] == 1:
						if strand == "+":
							lout.append("alt_3p")
						else:
							lout.append("alt_5p")
					else:
						if strand == "+":
							lout.append("alt_end")
						else:
							lout.append("alt_start")
	
					if right_alts[i]["novel"] == 1:
						lout[-1] += "*"
	
					lout.append(",".join(right_alts[i]["A_tid"]))
					lout.append(",".join(right_alts[i]["B_tid"]))
					
					lout.append(",".join(map(loc_reformat, right_alts[i]["A"])))
					lout.append(",".join(map(loc_reformat, right_alts[i]["B"])))
					
					print "\t".join(lout)
			
			for i in range(len(alt_cas)):
				
				print_me = True
				if args.no_novel_vs_novel:
					if alt_cas[i]["a_novel"] == True and alt_cas[i]["b_novel"] == True:
						print_me = False
				
				if print_me:
				
					aloc = "ALOC_{:08d}".format(loc_index)
					loc_index += 1
					
					lout = [aloc, gid, gname, strand, "alt_cassette:" + str(len(alt_cas[i]["A"])-1)]
	
					if alt_cas[i]["novel"] == 1:
						lout[-1] += "*"
					
					lout.append(",".join(alt_cas[i]["A_tid"]))
					lout.append(",".join(alt_cas[i]["B_tid"]))
					
					lout.append(",".join(map(loc_reformat, alt_cas[i]["A"])))
					lout.append(",".join(map(loc_reformat, alt_cas[i]["B"])))
					
					print "\t".join(lout)
			
			for i in range(len(me_pairs)):
				
				print_me = True
				if args.no_novel_vs_novel:
					if me_pairs[i]["a_novel"] == True and me_pairs[i]["b_novel"] == True:
						print_me = False
				
				if print_me:
				
					if me_pairs[i]["special"] == 0:
						aloc = "ALOC_{:08d}".format(loc_index)
						loc_index += 1
		
						lout = [aloc, gid, gname, strand, "mutually_exclusive"]
						if me_pairs[i]["novel"] == 1:
							lout[-1] += "*"
							
						lout.append(",".join(me_pairs[i]["A_tid"]))
						lout.append(",".join(me_pairs[i]["B_tid"]))
		
						lout.append(",".join(map(loc_reformat, me_pairs[i]["A"])))
						lout.append(",".join(map(loc_reformat, me_pairs[i]["B"])))
						
						print "\t".join(lout)

			for i in range(len(retained_introns)):

				print_me = True
				if args.no_novel_vs_novel:
					if retained_introns[i]["a_novel"] == True and retained_introns[i]["b_novel"] == True:
						print_me = False
				
				if print_me:
					
					aloc = "ALOC_{:08d}".format(loc_index)
					loc_index += 1
	
					lout = [aloc, gid, gname, strand, "retained_intron"]
					if retained_introns[i]["novel"] == 1:
						lout[-1] += "*"
						
					lout.append(",".join(retained_introns[i]["A_tid"]))
					lout.append(",".join(retained_introns[i]["B_tid"]))
	
					lout.append(",".join(map(loc_reformat, retained_introns[i]["A"])))
					lout.append(",".join(map(loc_reformat, retained_introns[i]["B"])))
					
					print "\t".join(lout)
												
	return 0

def find_retained_introns(introns, exons, itoe):
	# searches locus for retained introns..
	#
	# <----[===============]---->
	# <----[=]-----------[=]---->
	#
	# these will have to be quantifed as the junction verses 
	# reads that map over the junction's location so the
	# output "alt" regions need to be windows around the 
	# junction's two anchors.  We're looking for any 
	# exon that completely contains any intron. exons that
	# contain more than one intron will be considered more
	# than on retained intron event.
	#
	
	rwin = 8 # window size for checking for inclusion of a retained intron from alignments
	alts = []
	used_introns = set([])
	
	# check each exon against the intron table
	for eid in exons:
		efeature = eid.split("|")
		efi = [int(efeature[1]), int(efeature[2])]
		for iid in introns:
			if iid not in used_introns:
				ifeature = iid.split("|")
				ifi = [int(ifeature[1]), int(ifeature[2])]
				
				if efi[0] < ifi[0] and efi[1] > ifi[1]:
					# got one
					temp = new_alt()
					# the alt will be the junction
					temp["B"].append(iid)
					
					# gather the length of the intron and the lengths of the two
					# exons on either side of it to make sure the rwin window
					# is feasable. if not, adjust it.
					intron_length = feature_length_from_id(iid)
					exon_left_length = 0
					exon_right_length = 0
					iexons = list(itoe[iid])
					for eid in iexons:
						if get_edge(eid, 1) < get_edge(iid, 0):
							exon_left_length = feature_length_from_id(eid)
						else:
							exon_right_length = feature_length_from_id(eid)
					
					# two windows, a and b, with left and right lengths from center
					win_a_left = rwin
					win_a_right = rwin
					win_b_left = rwin
					win_b_right = rwin
					if exon_left_length < rwin:
						win_a_left = exon_left_length
					if intron_length < rwin:
						win_a_right = intron_length-1
						win_b_left = intron_length-1
					if exon_right_length < rwin:
						win_b_right = exon_right_length
					
					temp["A"].append("|".join(map(str, [ifeature[0], ifi[0]-win_a_left, ifi[0]+win_a_right])))
					temp["A"].append("|".join(map(str, [ifeature[0], ifi[1]-win_b_left, ifi[1]+win_b_right])))
					
					# copy transcripts for the intron and the exon
					temp["A_tid"] = list(exons[eid])
					temp["B_tid"] = list(introns[iid])
					
					# check if the transcript list for either A or B is completely novel
										
#					type = set([])
#					for tid in temp["A_tid"]:
#						if _NOVEL_PREFIX in tid:
#							type.update(["n"])
#						else:
#							type.update(["a"])
#					
#					if len(type) == 1 and "n" in type:
#						temp["novel"] = 1
#					else:
#						# check B form
#						type = set([])
#						for tid in temp["B_tid"]:
#							if _NOVEL_PREFIX in tid:
#								type.update(["n"])
#							else:
#								type.update(["a"])
#						
#						if len(type) == 1 and "n" in type:
#							temp["novel"] = 1 
#					
					alts.append(temp)
					used_introns.update([iid])

	# update the novel status for each alt location
	for i in range(len(alts)):

		temp_a = set([])
		for tid in alts[i]["A_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_a.update(["n"])
			else:
				temp_a.update(["a"])
		
		temp_b = set([])
		for tid in alts[i]["B_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_b.update(["n"])
			else:
				temp_b.update(["a"])
		
		alts[i]["a_novel"] = (len(temp_a) == 1 and "n" in temp_a)
		alts[i]["b_novel"] = (len(temp_b) == 1 and "n" in temp_b)
				
		# exclusive or to avoid novel on novel comparisons
		if (len(temp_a) == 1 and "n" in temp_a) or (len(temp_b) == 1 and "n" in temp_b):
			alts[i]["novel"] = 1
					
	return alts

def find_me_pairs(introns, exons, itoe, ti):
	# searches for mutually exclusive exon arrangments between two introns:
	# []-------------------[]--------[]
	# []-----------[]----------------[]
	#
	# see Evernote notes - this should be easy
	#
	# A = long,short
	# B = short, long
	#
	
	alts = []
	left_iid = {}
	left_tid = {}	
	right_iid = {}
	intron_index = 0
	made_hash = {}

	# put left anchors into buckets. store transcript id and intron id
	for iid in introns.keys():
		aid = get_lanchor(iid)

		for tid in introns[iid]:
			if len(ti[tid]) > 1:
				# only consider introns from transcripts that have at least 2 introns
				
				if aid not in left_iid:
					left_iid[aid] = []
					left_tid[aid] = []
					right_iid[aid] = []
				
				left_tid[aid].append(tid)
				left_iid[aid].append(iid)
				right_iid[aid].append("")
	
	# loop through the left anchors. if there is one with more than one unique transcript
	# id in the list then test if for any pair of those the next intron's right 
	# anchors line up between the two.
	
	for aid in left_tid.keys():
		if len(left_tid[aid]) == 1:
			# only one transcript at this anchor, pop it out
			left_tid.pop(aid)
			left_iid.pop(aid)
			right_iid.pop(aid)
			
		else:
			# more than one, check each transcript's intron to see if there is a next intron 
			# at all - if not remove them
			
			i = len(left_tid[aid]) - 1
			while i >= 0:
				tid = left_tid[aid][i]
				iid = left_iid[aid][i]
				intron_index = -1
				
				for j in range(len(ti[tid])):
					if ti[tid][j] == iid:
						intron_index = j
						break
				
				if j+1 == len(ti[tid]):
					# last intron remove this one from the set
					left_tid[aid].pop(i)
					left_iid[aid].pop(i)
					right_iid[aid].pop(i)
				
				else:
					# snag the right side intron
					right_iid[aid][i] = ti[tid][j+1]
				
				i -= 1
			
			# check remaining length
			if len(left_tid[aid]) < 2:
				# remove
				left_tid.pop(aid)
				left_iid.pop(aid)
				right_iid.pop(aid)
			
	# if anything is left can check for mutually exclusive splicing conditions
	if len(left_tid.keys()) > 0:
		
		for aid in left_tid.keys():
#			print aid
			# loop through these left anchors. for any pair check if the right anchors
			# of the NEXT introns match up. this forms a mutually exclusive splice
			# configuration.
			n = len(left_tid[aid])
			for i in range(n-1):
				for j in range(i+1, n):
					# the right anchors have to match, the right introns can't be identical and they can't
					# be from the same transcript
					if get_ranchor(right_iid[aid][i]) == get_ranchor(right_iid[aid][j]) and right_iid[aid][i] != right_iid[aid][j] and left_iid[aid][i] != left_iid[aid][j]:
						# got one

						temp = new_alt()
						if feature_length_from_id(left_iid[aid][i]) > feature_length_from_id(left_iid[aid][j]):
							# A = long, short
							# B = short, long
							temp["A"] = [left_iid[aid][i], right_iid[aid][i]]
							temp["B"] = [left_iid[aid][j], right_iid[aid][j]]
						else:
							temp["B"] = [left_iid[aid][i], right_iid[aid][i]]
							temp["A"] = [left_iid[aid][j], right_iid[aid][j]]
						
						hid = ",".join(temp["A"]) + ";" + ",".join(temp["B"])
						
						if hid not in made_hash:
							made_hash[hid] = 1
							
							# figure out if this is a *special* type where it looks like a 
							# mutually exclusive event but both exons are included together
							# in another isoform. 

							# get the eid's of the two exons then intersect their transcript lists.
							eid_a = list(itoe[temp["A"][0]].intersection(itoe[temp["A"][1]]))
							eid_b = list(itoe[temp["B"][0]].intersection(itoe[temp["B"][1]]))
							
							if len(eid_a) == 1 and len(eid_b) == 1:
								# continue
								seid_a = set(exons[eid_a[0]])
								seid_b = set(exons[eid_b[0]])
								tid_int = seid_a.intersection(seid_b)
								
								if len(tid_int) > 0:
									temp["special"] = 1
								
								alts.append(temp)
							else:
								sys.stderr.write("[find_me_pairs] found mutually exclusive pair somehow missing exons\n")

	
	# for each completed alt find the transcript lists for the ref and alt features
	for i in range(len(alts)):
		iid_a = alts[i]["A"]
		iid_b = alts[i]["B"]
		
		tid_a = set(list(introns[iid_a[0]])) & set(list(introns[iid_a[1]]))
		tid_b = set(list(introns[iid_b[0]])) & set(list(introns[iid_b[1]]))
		
		alts[i]["A_tid"] = list(tid_a)
		alts[i]["B_tid"] = list(tid_b)

		temp_a = set([])
		for tid in alts[i]["A_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_a.update(["n"])
			else:
				temp_a.update(["a"])
		
		temp_b = set([])
		for tid in alts[i]["B_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_b.update(["n"])
			else:
				temp_b.update(["a"])
		
		if (len(temp_a) == 1 and "n" in temp_a) or (len(temp_b) == 1 and "n" in temp_b):
			alts[i]["novel"] = 1

		alts[i]["a_novel"] = (len(temp_a) == 1 and "n" in temp_a)
		alts[i]["b_novel"] = (len(temp_b) == 1 and "n" in temp_b)


	return alts
	

def find_alt_cassettes(introns, ti):
	# looks through introns in the hash. for each intron lists of intron ids
	# paired with their transcript ids are made for each of the intron's two 
	# anchors including only introns that are shorter than it. if there are 
	# pairs of transcript ids between the two ends then we have a cassette. the
	# final cassette with have the ref intron and ALL introns between its 
	# anchors from the other transcript. if this happens with more than one
	# transcript then more than one pairing will be made
	# 
	# []-----------[]------------[]
	# []-------------------------[]
	#
	# or
	#
	# []-------[]--------[]------[]
	# []-------------------------[]
	#
	# or more!
	#
	
	alts = []
	iids = introns.keys()
	left_id = ""
	right_id = ""
	made_hash = {}

	if len(iids) < 2:
		return alts
	else:
		# continue
		
		for iid in iids:
			left_iid = []
			left_tid = []
			right_iid = []
			right_tid = []
			tid_int = set([])
			
			left_id = get_lanchor(iid)
			right_id = get_ranchor(iid)
			ilen = feature_length_from_id(iid)
			
			for i in range(len(iids)):
				if iids[i] != iid and ilen > feature_length_from_id(iids[i]):
					# intron is a potential, check anchors
					if left_id == get_lanchor(iids[i]):
						# append the intron id and its transcript information
						for tid in introns[iids[i]]:
							left_iid.append(iids[i])
							# tids are concatenated with the intron's index, split that out
							# left_tid.append(tid.split(":")[0])
							left_tid.append(tid)

					if right_id == get_ranchor(iids[i]):
						# append the intron id and its transcript information
						for tid in introns[iids[i]]:
							right_iid.append(iids[i])
							# tids are concatenated with the intron's index, split that out
							# right_tid.append(tid.split(":")[0])
							right_tid.append(tid)
			
			
			# intersect the left and right tid lists
			tid_int = set(left_tid).intersection(set(right_tid))

			# if there is an intersection we have more to do
			if len(tid_int) > 0:
				
				for tid in list(tid_int):
					# track down the intron ids that go with this transcript id
					lid = ""
					rid = ""
					for i in range(len(left_tid)):
						if left_tid[i] == tid:
							lid = left_iid[i]
							break
					for i in range(len(right_tid)):
						if right_tid[i] == tid:
							rid = right_iid[i]
							break
					
					# now we almost have all we need
					#alts.append(new_alt())
					#alts[-1]["ref"] = iid
					temp = new_alt()
					temp["B"].append(iid)
#					print lid, rid
					
					# loop through the transcript's intron list and append each one
					# between the two found above, including the two found. the
					# introns should be sorted so this should work.
					
					i = 0
					keep = False
					while i < len(ti[tid]):
						if ti[tid][i] == lid:
							keep = True
#							print "started"
							temp["A"].append(ti[tid][i])
						elif ti[tid][i] == rid:
#							print "found right id"
							temp["A"].append(ti[tid][i])
							break
						elif keep:
							temp["A"].append(ti[tid][i])

#						print i
						i += 1
					
					# check if this alt cassette is already in the table
					hkey = temp["B"][0] + ";" + ",".join(temp["A"])
					if hkey not in made_hash:
						made_hash[hkey] = 1
						alts.append(temp)

	# for each completed alt find the transcript lists for the ref and alt features
	for i in range(len(alts)):
		iid_a = alts[i]["A"]
		iid_b = alts[i]["B"]
		
		tid_a = set(list(introns[iid_a[0]]))		
		tid_b = set(list(introns[iid_b[0]]))
		
		if len(iid_a) > 1:
			for j in range(1, len(iid_a)):
				tid_a &= set(list(introns[iid_a[j]]))
				
		alts[i]["B_tid"] = list(tid_b)
		alts[i]["A_tid"] = list(tid_a)
		
		# check for novel
		temp_a = set([])
		for tid in alts[i]["A_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_a.update(["n"])
			else:
				temp_a.update(["a"])
		
		temp_b = set([])
		for tid in alts[i]["B_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_b.update(["n"])
			else:
				temp_b.update(["a"])
		
		# exclusive or to avoid novel on novel tests
		if (len(temp_a) == 1 and "n" in temp_a) or (len(temp_b) == 1 and "n" in temp_b):
			alts[i]["novel"] = 1

		alts[i]["a_novel"] = (len(temp_a) == 1 and "n" in temp_a)
		alts[i]["b_novel"] = (len(temp_b) == 1 and "n" in temp_b)

	return alts

# --
# find_alt_ends
# from a list of "end" introns, either left or right side, and a full list of 
# all introns, this function finds pairings between alternative ends.
#
def find_alt_ends(end_introns, introns, end_exons, etoi, itoe, mode):
	# finds these...
	#
	# []----------------------------[]------->
	#          []-------------------[]------->
	# or
	# <----[]----------------------------[]------->
	#               []-------------------[]------->
	#
	# The rule is that the "end" intron, from end_introns, must share 
	# the first anchor with another transcript and the start exon
	# must not be shared. since the tables only have unique introns
	# there is no check for left and right anchor alignment between
	# two introns.
	#
	# Also need this category:
	#
	#              a
	# [   ]---------------------[]---->
	# [  ]----------------------[]---->
	#              b
	#
	# For left side this is equal right anchors and unequal left anchors but 
	# the left anchors hit overlapping exons. Measure the distance of the
	# overlap of b past the left anchor of a. Report 'a' as ref and 'b' as the alt.
	#
	
	alts = []
	aset = {}
	aprimes = {}
	intron_exons = {}
	iidlist = end_introns.keys()
	alt_primes = set([])
	
	if mode == "left":
		# left introns compare right anchors
		end_id = 2
	else:
		# right introns compare left anchors
		end_id = 1
	
	if len(introns) < 2:
		return alts

	# make a full set of anchors from the end introns
	for iid in end_introns:
		intron = iid.split("|")
		aid = intron[end_id]
		
		if aid not in aset:
			aset[aid] = []
		
		aset[aid].append(iid)			
	
	# search out pairs of introns that anchor to overlapping exons that do not share
	# the same edge (aka the alt-3' and alt-5' types rather than the alt-start and 
	# alt-end types where the end exons do not overlap at all)

	for aid in aset.keys():
		
		if len(aset[aid]) > 1:
			iidlist = aset[aid]
			n = len(iidlist)
			alt_primes = set([])
			
			# pull together the end exons associated with each intron
			# in this bucket
			for iid in aset[aid]:
				intron_exons[iid] = []
				elist = list(itoe[iid])
				for eid in elist:
					if eid in end_exons:
						intron_exons[iid].append(eid)
					
			# look for pairings of introns that have overlapping exons at their left or right end
			for i in range(n-1):
				elist1 = []
				for eid in intron_exons[iidlist[i]]:
					elist1.append(eid.split("|"))
				
				if mode == "left":
					edge1 = elist1[0][2]
				else:
					edge1 = elist1[0][1]
				
				for j in range(i+1, n):
					elist2 = []
					for eid in intron_exons[iidlist[j]]:
						elist2.append(eid.split("|"))
		
					if mode == "left":
						edge2 = elist2[0][2]
					else:
						edge2 = elist2[0][1]
					
					# look for any feature overlap between elist1 and elist2 with
					# mismatched edges
					ovr = False
					if edge1 != edge2:
						h = 0
						while h < len(elist1) and not ovr:
							f1 = elist1[h]
							for k in range(len(elist2)):
								f2 = elist2[k]
								
								if feature_overlap(f1, f2):
									ovr = True
									break
							
							h += 1
					
					if ovr:
						# the introns at i and j both land into exons that overlap and do not
						# share edges on the intron side so this is probably an alt prime.
						# i and j should be removed from this aid bucket
						alt_primes.update([i, j])
						
						# we're gonna move introns i and j over to the aprimes table
						# and out of the aset table which will remain intact for 
						# finding alt start/ends
						
						if aid not in aprimes:
							aprimes[aid] = set([])
						
						aprimes[aid].update([iidlist[i], iidlist[j]])
			
			# NOTE: removing this next section because it seems that the analysis doesn't work 
			# correctly. these introns need to be compared for alt-prime status where they
			# are pitted in a 1 vs 1 type setup but they also need to be evaluated as alternate
			# start/end sites.
							
			# finished processing this aid, now clean out the introns that are moved
			# to the aprimes table
			#rem_idx_list = sorted(list(alt_primes))
			#n = len(rem_idx_list)
			#while n:
			#	n -= 1
			#	aset[aid].pop(rem_idx_list[n])
			
			
	# add in additional introns that might fall into the aset buckets. since introns
	# are unique in these tables there would not be a case of an end intron
	# sharing a right and left anchor with another intron - so there's no point
	# in checking for that condition.
	for iid in introns:
		intron = iid.split("|")
		aid = intron[end_id]
		
		if aid in aset and iid not in end_introns:
			aset[aid].append(iid)
			
	
	# now we have a list of end_introns sorted by their right anchors.
	# for any of them with 2 or more end_introns we can setup pairings
	for aid in aset:
		if len(aset[aid]) > 1:
			# more than one left side intron sharing right-side 
			# anchors.
			
			if len(aset[aid]) == 2:
				# only one pairing necessary make sure at least one of these two is in the
				# end_introns table. 
				if aset[aid][0] in end_introns or aset[aid][1] in end_introns: 
					alts.append(new_alt())
					# make sure that the intron that's the "end" intron is the ref intron
					if aset[aid][0] in end_introns:					
						alts[-1]["A"].append(aset[aid][0])
						alts[-1]["B"].append(aset[aid][1])
					else:
						alts[-1]["A"].append(aset[aid][1])
						alts[-1]["B"].append(aset[aid][0])
			else:
				# more than two so we'll setup pairings of each one verses each
				# of the rest for any possible ref's that are "end introns"
				for i in range(len(aset[aid])):
					if aset[aid][i] in end_introns:
						alts.append(new_alt())
						for j in range(len(aset[aid])):
							if i == j:
								alts[-1]["A"].append(aset[aid][j])
							else:
								alts[-1]["B"].append(aset[aid][j])

	# process the alt-primes
	for aid in aprimes:
		aprimes[aid] = list(aprimes[aid])
		if len(aprimes[aid]) > 1:
			# more than one left side intron sharing right-side 
			# anchors.
			
			if len(aprimes[aid]) == 2:
				# only one pairing necessary make sure at least one of these two is in the
				# end_introns table. 
				if aprimes[aid][0] in end_introns or aprimes[aid][1] in end_introns: 
					alts.append(new_alt())
					alts[-1]["special"] = 1
					
					# make the shorter intron the A form
					if feature_length_from_id(aprimes[aid][0]) < feature_length_from_id(aprimes[aid][1]):					
						alts[-1]["A"].append(aprimes[aid][0])
						alts[-1]["B"].append(aprimes[aid][1])
					else:
						alts[-1]["A"].append(aprimes[aid][1])
						alts[-1]["B"].append(aprimes[aid][0])
						
			else:
				# more than two so we'll setup pairings of each one verses each
				# of the rest for any possible ref's that are "end introns"
				for i in range(len(aprimes[aid])):
					alts.append(new_alt())
					alts[-1]["special"] = 1
					for j in range(len(aprimes[aid])):
						if i == j:
							alts[-1]["A"].append(aprimes[aid][j])
						else:
							alts[-1]["B"].append(aprimes[aid][j])

									
	# for each completed alt find the transcript lists for the ref and alt features
	for i in range(len(alts)):
		iid_a = alts[i]["A"]
		iid_b = alts[i]["B"]
		
		# A always only has a single intron while B may have more than one
		
		tid_a = set(list(introns[iid_a[0]]))
		tid_b = set(list(introns[iid_b[0]]))
		
		if len(iid_b) > 1:
			for j in range(1, len(iid_b)):
				tid_b |= set(list(introns[iid_b[j]]))
		
		alts[i]["A_tid"] = list(tid_a)
		alts[i]["B_tid"] = list(tid_b)

		temp_a = set([])
		for tid in alts[i]["A_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_a.update(["n"])
			else:
				temp_a.update(["a"])
		
		temp_b = set([])
		for tid in alts[i]["B_tid"]:
			if _NOVEL_PREFIX in tid:
				temp_b.update(["n"])
			else:
				temp_b.update(["a"])
		
		if (len(temp_a) == 1 and "n" in temp_a) or (len(temp_b) == 1 and "n" in temp_b):
			alts[i]["novel"] = 1

		alts[i]["a_novel"] = (len(temp_a) == 1 and "n" in temp_a)
		alts[i]["b_novel"] = (len(temp_b) == 1 and "n" in temp_b)

	# return result
	return alts

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
parser.add_argument('--no-novel-vs-novel', dest="no_novel_vs_novel", action="store_const", const=True, default=False, 
				help="Disallow novel vs novel alternative splice locations (default: allowed)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
