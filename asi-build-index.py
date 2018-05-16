#!/usr/bin/env python
#==============================================================================
# main.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Working on a graph based alt-splice indexing thing 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser, basename
import igraph as ig

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

# gtf fields
_CHROM = 0
_STRAND = 6
_START = 3
_END = 4
_TYPE = 2
_ATTR = 8
STRAND_POS = 1
STRAND_NEG = 0

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	dchroms = {}
	
	# check inputs
	if not isfile(args.gtf):
		sys.stderr.write("Error: input GTF does not exist ({})".format(args.gtf))
		return 1

	# build output file name
	fout_name = re.sub("\.gtf$", ".asrefbeta", basename(args.gtf))

	# load the gtf
	dchroms = load_gtf(args.gtf)

	gnodes = {}
	eset = {}
	iset = {}
	gexons = []
	gintrons = []
	# list of edge descriptions which will include the transcript names each edge belongs to, 
	# if it is an exon or not and its start/end coordinates
	gedges = []
	# list of edges in terms of vertices - this is the list that will be added to the graph
	ledge_info = []
	num_verts = 0
	dtid_nodes = {}
	aloc_index = 0
	pid_index = 0

	# open output file
	fout = open(fout_name, "w")

	g = None

	for k in sorted(dchroms.keys()):
		
		# per chrom we have a mess of transcripts. we will build the gnodes, eset and iset
		# dicts from them now
		gnodes = {}
		eset = {}
		iset = {}
		gedges = []
		ledge_info = []
		dtid_nodes = {}
		tid_strand = {}
		tid_gid = {}
		tid_gname = {}

		for tid in dchroms[k].keys():
			transcript_to_ve(dchroms[k][tid], gnodes, eset, iset)
			tid_strand[tid] = get_strand(dchroms[k][tid][0])
			tid_gid[tid] = get_gid(dchroms[k][tid][0])
			tid_gname[tid] = get_gname(dchroms[k][tid][0])

		# done with that, now eset and iset must be converted the main edge list
		# for the graph

		for eid in eset:
			ltmp = eid.split(":")			
			gedges.append((gnodes[ltmp[0]][0], gnodes[ltmp[1]][0]))
			# make sure the info has the correct genome oriented order of the two 
			# end point nodes
			if int(ltmp[0]) > int(ltmp[1]):
				ledge_info.append([ltmp[1], ltmp[0], "exon", eset[eid], STRAND_NEG])
			else:
				# if the left and right coordinates are already in the correct order then this was
				# probably from a positive strand feature
				ledge_info.append([ltmp[0], ltmp[1], "exon", eset[eid], STRAND_POS])

		for iid in iset:
			ltmp = iid.split(":")
			gedges.append((gnodes[ltmp[0]][0], gnodes[ltmp[1]][0]))
			# make sure the info has the correct genome oriented order of the two 
			# end point nodes
			if int(ltmp[0]) > int(ltmp[1]):
				ledge_info.append([ltmp[1], ltmp[0], "intron", iset[iid], STRAND_NEG])
			else:
				ledge_info.append([ltmp[0], ltmp[1], "intron", iset[iid], STRAND_POS])

		# get number of vertices
		num_verts = len(gnodes.keys())

		lverts = ["" for i in range(num_verts)]
		for vid in gnodes.keys():
			lverts[gnodes[vid][0]] = vid

		# make graph
		g = ig.Graph(directed=True)
		g.add_vertices(num_verts)
		g.add_edges(gedges)

		# annotate the vertices:
		for vid in gnodes.keys():
			g.vs[gnodes[vid][0]]["name"] = "{}:{}".format(k, vid)
			g.vs[gnodes[vid][0]]["tid_set"] = list(gnodes[vid][1])

		# annotate the edges
		for i in range(len(ledge_info)):
			g.es[i]['name'] = "{}:{}".format(k, "-".join(ledge_info[i][0:2]))
			g.es[i]['start'] = ledge_info[i][0]
			g.es[i]['end'] = ledge_info[i][1]
			g.es[i]['type'] = ledge_info[i][2]
			g.es[i]['tid_set'] = ledge_info[i][3]
			g.es[i]['strand'] = ledge_info[i][4]

		# now build a dict of transcript ids containing the vertices for each
		for i in range(num_verts):
			for tid in g.vs[i]['tid_set']:
				if tid not in dtid_nodes:
					dtid_nodes[tid] = []

				# append this node index to the transcript's dict
				dtid_nodes[tid].append(g.vs[i]['name'])

		# cluster it
		gi = g.clusters(mode=ig.WEAK)

		# groups
		groups = gi.membership

		# make a dict with cluster memberships
		dgroups = {}
		for i in range(len(groups)):
			if groups[i] not in dgroups:
				dgroups[groups[i]] = []

			dgroups[groups[i]].append(i)

		for gid in sorted(dgroups.keys()):

			# each cluster may contain alternative events
			ghat = g.subgraph(dgroups[gid])

			# --
			# 
			# Run alt cassette search
			#
			# --

			lres = find_alt_cassettes(ghat, dtid_nodes)

			if len(lres) > 0:
				# output alternative events
				for i in range(len(lres)):
					# iterate over the alternative path bundles
					aid = "ALOC_{:08d}".format(aloc_index)
					aloc_index += 1
					for j in range(len(lres[i])):
						# iterate to each path in the bundle
						pid = "PID_{:08d}".format(pid_index)
						pid_index += 1
						edge_list = lres[i][j][0].split("|")
						tid_list = lres[i][j][1]
						strand = tid_strand[tid_list[0]]

						# output the aloc id, the path id, type, number of exons, the intron set
						intron_list = []
						k = 0
						while k < len(edge_list):
							intron_list.append(trim_feature(edge_list[k]))
							k += 2

						lout = [aid, 
							pid, 
							tid_gid[tid_list[0]], 
							tid_gname[tid_list[0]], 
							"alt.cassette", 
							"{}".format((len(edge_list)-1)/2), 
							strand, 
							"|".join(intron_list), 
							",".join(lres[i][j][1])]

						fout.write("\t".join(lout))
						fout.write("\n")

			# --
			# 
			# Run alt-3' search
			#
			# --

			lres = find_alt_3p(ghat, dtid_nodes)

			if len(lres) > 0:
				# output alternative events
				for i in range(len(lres)):
					# iterate over the alternative path bundles
					aid = "ALOC_{:08d}".format(aloc_index)
					aloc_index += 1
					for j in range(len(lres[i])):
						# iterate to each path in the bundle
						pid = "PID_{:08d}".format(pid_index)
						pid_index += 1
						tid_list = lres[i][j][1]
						strand = tid_strand[tid_list[0]]

						lout = [aid, 
							pid, 
							tid_gid[tid_list[0]], 
							tid_gname[tid_list[0]], 							
							"alt.3p", 
							"0",
							strand, 
							trim_feature(lres[i][j][0]), 
							",".join(lres[i][j][1])]

						fout.write("\t".join(lout))
						fout.write("\n")

			# --
			# 
			# Run alt-5' search
			#
			# --

			lres = find_alt_5p(ghat, dtid_nodes)

			if len(lres) > 0:
				# output alternative events
				for i in range(len(lres)):
					# iterate over the alternative path bundles
					aid = "ALOC_{:08d}".format(aloc_index)
					aloc_index += 1
					for j in range(len(lres[i])):
						# iterate to each path in the bundle
						pid = "PID_{:08d}".format(pid_index)
						pid_index += 1

						tid_list = lres[i][j][1]
						strand = tid_strand[tid_list[0]]

						lout = [aid, 
							pid, 
							tid_gid[tid_list[0]], 
							tid_gname[tid_list[0]], 							
							"alt.5p", 
							"0",
							strand, 
							trim_feature(lres[i][j][0]), 
							",".join(lres[i][j][1])]

						fout.write("\t".join(lout))
						fout.write("\n")

			# --
			# 
			# find mutually exclusive 
			#
			# --

			lres = find_mut_ex(ghat, dtid_nodes)

			if len(lres) > 0:
				# output alternative events
				for i in range(len(lres)):
					# iterate over the alternative path bundles
					aid = "ALOC_{:08d}".format(aloc_index)
					aloc_index += 1
					for j in range(len(lres[i])):
						# iterate to each path in the bundle
						pid = "PID_{:08d}".format(pid_index)
						pid_index += 1

						tid_list = lres[i][j][1]
						strand = tid_strand[tid_list[0]]
						edge_list = lres[i][j][0].split("|")
						for k in range(len(edge_list)):
							edge_list[k] = trim_feature(edge_list[k])						

						lout = [aid, 
							pid, 
							tid_gid[tid_list[0]], 
							tid_gname[tid_list[0]], 							
							"mut.ex", 
							"0",
							strand, 
							"|".join(edge_list), 
							",".join(lres[i][j][1])]

						fout.write("\t".join(lout))
						fout.write("\n")
						#print lout

			# --
			# 
			# run alt start search
			#
			# --

			lres = find_alt_starts(ghat, dtid_nodes)

			if len(lres) > 0:
				# output alternative events
				for i in range(len(lres)):
					# iterate over the alternative path bundles
					aid = "ALOC_{:08d}".format(aloc_index)
					aloc_index += 1
					for j in range(len(lres[i])):
						# iterate to each path in the bundle
						pid = "PID_{:08d}".format(pid_index)
						pid_index += 1

						tid_list = lres[i][j][1]
						strand = tid_strand[tid_list[0]]

						lout = [aid, 
							pid, 
							tid_gid[tid_list[0]], 
							tid_gname[tid_list[0]], 							
							"alt.start", 
							"0",
							strand, 
							trim_feature(lres[i][j][0]), 
							",".join(lres[i][j][1])]

						fout.write("\t".join(lout))
						fout.write("\n")

			# --
			# 
			# run alt end search
			#
			# --

			lres = find_alt_ends(ghat, dtid_nodes)

			if len(lres) > 0:
				# output alternative events
				for i in range(len(lres)):
					# iterate over the alternative path bundles
					aid = "ALOC_{:08d}".format(aloc_index)
					aloc_index += 1
					for j in range(len(lres[i])):
						# iterate to each path in the bundle
						pid = "PID_{:08d}".format(pid_index)
						pid_index += 1

						tid_list = lres[i][j][1]
						strand = tid_strand[tid_list[0]]

						lout = [aid, 
							pid, 
							tid_gid[tid_list[0]], 
							tid_gname[tid_list[0]], 							
							"alt.end", 
							"0",
							strand, 
							trim_feature(lres[i][j][0]), 
							",".join(lres[i][j][1])]

						fout.write("\t".join(lout))
						fout.write("\n")


	fout.close()

	return 0


# --
# find_alt_cassettes
# see NOTES.txt. this function will search out all alternative cassette type
# path bundles
def find_alt_cassettes(g, dtid):

	lsets = []
	all_introns = False
	dpaths = {}


	# get all output lists - first is the list, per node, of adjacent nodes
	# and the second is the same but the inner list elements are edges
	lvout = g.get_adjlist(mode=ig.OUT)
	lvin = g.get_adjlist(mode=ig.IN)
	leout = g.get_inclist(mode=ig.OUT)
	lein = g.get_inclist(mode=ig.IN)


	for i in range(len(lvout)):
		if len(lvout[i]) > 1:
			# this node has multiple outputs - are they introns?
			all_introns = True
			for j in leout[i]:
				if g.es[j]['type'] != "intron":
					all_introns = False
					break

			if all_introns:
				# for any that connect to a node with alternative inputs...
				for j in range(len(lvout[i])):
					vid = lvout[i][j]
					eid = leout[i][j]
					# does this target node have multiple inputs?
					if len(lvin[vid]) > 1:
						# get tid sets for these two and the intron
						don_tid = set(g.vs[i]['tid_set'])
						acc_tid = set(g.vs[vid]['tid_set'])
						intron_tid = set(g.es[eid]['tid_set'])
						don_name = g.vs[i]['name']
						acc_name = g.vs[vid]['name']

						don_tid = don_tid - intron_tid
						acc_tid = acc_tid - intron_tid
						shared_tid = don_tid & acc_tid

						if len(shared_tid) > 0:
							# we have at least one transcript that shares this 
							# donor and acceptor. now get the paths from those
							# from the donor to the acceptor

							shared_tid = list(shared_tid)
							dpaths = {}
							for tid in shared_tid:
								# get this transcript's subgraph
								vid_set = dtid[tid]
								g_tid = g.subgraph(vid_set)
								# remove any edges that do not belong to this tid
								edrops = []
								for k in range(len(g_tid.es)):
									if tid not in g_tid.es[k]['tid_set']:
										edrops.append(k)

								# drop edges
								if len(edrops) > 0:
									g_tid.delete_edges(edrops)

								ltmp = g_tid.get_shortest_paths(don_name, to=acc_name, mode=ig.OUT, output="epath")
								# there should only be one path because this graph is only one transcript
								ltmpHat = [g_tid.es[k]['name'] for k in ltmp[0]]
								pid = "|".join(ltmpHat)
								if pid not in dpaths:
									dpaths[pid] = []
								dpaths[pid].append(tid)

							if len(dpaths.keys()) > 0:
								# we have at least one alternative path, build up the alt set which 
								# will be a list of each path including the initial single-intron path
								ltmp = [[g.es[eid]['name'], g.es[eid]['tid_set']]]
								for pid in dpaths.keys():
									ltmp.append([pid, dpaths[pid]])

								lsets.append(ltmp)


#	for i in range(len(lsets)):
#		print lsets[i]
						
	return lsets


# --
# find_alt_3p
# see NOTES.txt
# find all alternative 3' features (alt donor site with common acceptor)
def find_alt_3p(g, dtid):

	lsets = []
	all_exons = False
	dpaths = {}

	# get the chrom name
	ntmp = g.vs[0]['name']
	chrom = ntmp.split(":")[0]

	# get all output lists - first is the list, per node, of adjacent nodes
	# and the second is the same but the inner list elements are edges

	lvout = g.get_adjlist(mode=ig.OUT)
	lvin = g.get_adjlist(mode=ig.IN)
	leout = g.get_inclist(mode=ig.OUT)
	lein = g.get_inclist(mode=ig.IN)

	for i in range(len(lvout)):
		if len(lvout[i]) > 1:
			# this node has multiple outputs - are they all exons?
			all_exons = True
			for j in leout[i]:
				if g.es[j]['type'] != "exon":
					all_exons = False
					

			if not all_exons:
				continue

			# tid set
			ltid = g.vs[i]['tid_set']
			name0 = g.vs[i]['name']

			# build a dict keyed by the paths found
			dpaths = {}
			# per transcript id find the two step path by vertex names
			for tid in ltid:
				sz_path = name0
				vid_set = dtid[tid]
				g_tid = g.subgraph(vid_set)
				# remove any edges that do not belong to this tid
				edrops = []
				for k in range(len(g_tid.es)):
					if tid not in g_tid.es[k]['tid_set']:
						edrops.append(k)

				# drop edges
				if len(edrops) > 0:
					g_tid.delete_edges(edrops)

				# find the first node
				lvidx = [g_tid.vs.find(name=name0).index]
				g_tid_lout = g_tid.get_adjlist(mode=ig.OUT)
				
				for j in range(2):
					# get next and next
					if len(g_tid_lout[lvidx[j]]) > 0:
						lvidx.append(g_tid_lout[lvidx[j]][0])
					else:
						break

				if len(lvidx)==3:
					pid = "|".join(g_tid.vs[lvidx]['name'])
					if pid not in dpaths:
						dpaths[pid] = []

					dpaths[pid].append(tid)

			lpid = dpaths.keys()
			if len(lpid) > 1:

				# there is more than one path out of the initial node. now we have to figure out
				# if any of the paths for the alt condition we are looking for
				# use a dict to bundle 
				d = {}
				for j in range(len(lpid)):
					ltmp = lpid[j].split("|")

					if ltmp[2] not in d:
						d[ltmp[2]] = []

					d[ltmp[2]].append(j)

				# at this point for any key in d, if the length of the list of indices is 
				# greater than 1 then that's an alt set

				for vid in d.keys():
					if len(d[vid]) > 1:
						# this is a set. append the intron and the transcript set of each
						# to the lsets list. first make the intron
						ltmp_sets = []
						for j in range(len(d[vid])):
							ltmp = lpid[d[vid][j]].split("|")
							# intron of interest is between ltmp[1] and ltmp[2] but due to 
							# strand we might not know which order
							ltmpA = ltmp[1].split(":")
							ltmpB = ltmp[2].split(":")
							
							if int(ltmpA[1]) > int(ltmpB[1]):
								eid = "{}-{}".format(ltmp[2], ltmpA[1])
							else:
								eid = "{}-{}".format(ltmp[1], ltmpB[1])

							ltmp_sets.append([eid, dpaths[lpid[j]]])

						if len(ltmp_sets) > 0:
							lsets.append(ltmp_sets)
						
	return lsets


# --
# find_alt_5p
# see NOTES.txt
# find all alternative 5' features 
def find_alt_5p(g, dtid):

	lsets = []
	all_exons = False
	dpaths = {}

	# get the chrom name
	ntmp = g.vs[0]['name']
	chrom = ntmp.split(":")[0]

	# get all output lists - first is the list, per node, of adjacent nodes
	# and the second is the same but the inner list elements are edges

	lvout = g.get_adjlist(mode=ig.OUT)
	lvin = g.get_adjlist(mode=ig.IN)
	leout = g.get_inclist(mode=ig.OUT)
	lein = g.get_inclist(mode=ig.IN)

	# we have to look for nodes with alternative inputs where those input edges are all 
	# exons.
	for i in range(len(lvin)):
		if len(lvin[i]) > 1:
			# this node has multiple outputs - are they all exons?
			all_exons = True
			for j in lein[i]:
				if g.es[j]['type'] != "exon":
					all_exons = False
					

			if not all_exons:
				continue

			# get the transcript set that shares the alternative node then visit 
			# each transcript to get the 2-step path leading to it.

			# tid set
			ltid = g.vs[i]['tid_set']
			name0 = g.vs[i]['name']

			# build a dict keyed by the paths found
			dpaths = {}
			# per transcript id find the two step path by vertex names
			for tid in ltid:
				# subgraph out this transcript
				vid_set = dtid[tid]
				g_tid = g.subgraph(vid_set)
				# remove any edges that do not belong to this transcript
				edrops = []
				for k in range(len(g_tid.es)):
					if tid not in g_tid.es[k]['tid_set']:
						edrops.append(k)
				# drop edges if necessary
				if len(edrops) > 0:
					g_tid.delete_edges(edrops)

				# find the first node
				lvidx = [g_tid.vs.find(name=name0).index]
				# get the input lists which will cause us to traverse backwards
				g_tid_lin = g_tid.get_adjlist(mode=ig.IN)
				
				for j in range(2):
					# get next and next
					if len(g_tid_lin[lvidx[j]]) > 0:
						lvidx.append(g_tid_lin[lvidx[j]][0])
					else:
						break

				if len(lvidx)==3:
					pid = "|".join(g_tid.vs[lvidx]['name'])
					if pid not in dpaths:
						dpaths[pid] = []

					dpaths[pid].append(tid)

			lpid = dpaths.keys()
			if len(lpid) > 1:

				# there is more than one path out of the initial node. now we have to figure out
				# if any of the paths for the alt condition we are looking for
				# use a dict to bundle 
				d = {}
				for j in range(len(lpid)):

					ltmp = lpid[j].split("|")

					if ltmp[2] not in d:
						d[ltmp[2]] = []

					d[ltmp[2]].append(j)

				# at this point for any key in d, if the length of the list of indices is 
				# greater than 1 then that's an alt set

				for vid in d.keys():
					if len(d[vid]) > 1:
						# this is a set. append the intron and the transcript set of each
						# to the lsets list. first make the intron
						ltmp_sets = []
						for j in range(len(d[vid])):
							ltmp = lpid[d[vid][j]].split("|")
							# intron of interest is between ltmp[1] and ltmp[2] but due to 
							# strand we might not know which order
							ltmpA = ltmp[1].split(":")
							ltmpB = ltmp[2].split(":")
							
							if int(ltmpA[1]) > int(ltmpB[1]):
								eid = "{}-{}".format(ltmp[2], ltmpA[1])
							else:
								eid = "{}-{}".format(ltmp[1], ltmpB[1])

							ltmp_sets.append([eid, dpaths[lpid[j]]])

						if len(ltmp_sets) > 0:
							lsets.append(ltmp_sets)
						
	return lsets


# --
# find_mut_ex
# see NOTES.txt
# find mutually exclusive paths
def find_mut_ex(g, dtid):

	lsets = []
	all_introns = False
	dpaths = {}

	# get the chrom name
	ntmp = g.vs[0]['name']
	chrom = ntmp.split(":")[0]

	# get all output lists - first is the list, per node, of adjacent nodes
	# and the second is the same but the inner list elements are edges

	lvout = g.get_adjlist(mode=ig.OUT)
	lvin = g.get_adjlist(mode=ig.IN)
	leout = g.get_inclist(mode=ig.OUT)
	lein = g.get_inclist(mode=ig.IN)

	for i in range(len(lvout)):
		if len(lvout[i]) > 1:
			# this node has multiple outputs - are they all introns?
			all_introns = True
			for j in leout[i]:
				if g.es[j]['type'] != "intron":
					all_introns = False
					
			if not all_introns:
				continue

			# tid set
			ltid = g.vs[i]['tid_set']
			name0 = g.vs[i]['name']

			# build a dict keyed by the paths found
			dpaths = {}
			# per transcript id find the three step path by vertex names
			for tid in ltid:
				vid_set = dtid[tid]
				g_tid = g.subgraph(vid_set)
				# remove any edges that do not belong to this tid
				edrops = []
				for k in range(len(g_tid.es)):
					if tid not in g_tid.es[k]['tid_set']:
						edrops.append(k)

				# drop edges
				if len(edrops) > 0:
					g_tid.delete_edges(edrops)

				# find the first node
				lvidx = [g_tid.vs.find(name=name0).index]
				g_tid_lout = g_tid.get_adjlist(mode=ig.OUT)
				for j in range(3):
					# get next and next
					if len(g_tid_lout[lvidx[j]]) > 0:
						lvidx.append(g_tid_lout[lvidx[j]][0])
					else:
						break

				if len(lvidx)==4:
					pid = "|".join(g_tid.vs[lvidx]['name'])
					if pid not in dpaths:
						dpaths[pid] = []

					dpaths[pid].append(tid)

			lpid = dpaths.keys()

			if len(lpid) > 1:

				#
				# the alt donor site has multiple paths of the correct length.
				# first we have to find if any have a common end node.
				# for those that have common end nodes we have to figure out if 
				# the exon contained by each is mutually exclusive which means they
				# must not overlap and they must not both appear in the same transcript
				#

				d = {}
				for j in range(len(lpid)):
					ltmp = lpid[j].split("|")
					vid = ltmp[3]
					if vid not in d:
						d[vid] = []

					d[vid].append(j)

				for vid in d.keys():
					if len(d[vid]) > 1:
						# this end node appears in at least 2 paths...
						# make a list of the center exons

						eset = []
						for j in d[vid]:
							ltmp = lpid[j].split("|")
							# exon is ltmp[1] and ltmp[2]
							ltmpA = ltmp[1].split(":")
							ltmpB = ltmp[2].split(":")
							if int(ltmpA[1]) > int(ltmpB[1]):
								# A is downstream of B
								eid = "{}-{}".format(ltmp[2], ltmpA[1])
							else:
								eid = "{}-{}".format(ltmp[1], ltmpB[1])
							eset.append(eid)

						# look at this in terms of pairs of the eset
						
						for j in range(len(eset)-1):
							eset_hat = [eset[j]]
							for k in range(j+1, len(eset)):
								ltmp_sets = []
								eset_hat.append(eset[k])

								# make a dict of the transcript ids associated with the exons
								eset_tset = {}
								for eid in eset_hat:
									eidx = g.es.find(name=eid).index
									tset = g.es[eidx]['tid_set']
									for tid in tset:
										if tid not in eset_tset:
											eset_tset[tid] = 0
										eset_tset[tid] += 1

								# check the transcript counts. if any transcript count is 
								# 2 then this is not a valid pair
								valid_pair = True
								for tid in eset_tset.keys():
									if eset_tset[tid]==2:
										valid_pair = False

								# if this is a valid pair then the exons must also not
								# overlap each other
								if valid_pair and not feature_overlap(eset[j], eset[k]):
									#print eset[j], eset[k]
									#print lpid[d[vid][j]], lpid[d[vid][k]]
									
									for n in [j, k]:

										# make introns for position j
										ltmp = lpid[d[vid][n]].split("|")
										iset = [combine_vid_to_feature(ltmp[0], ltmp[1]), 
											combine_vid_to_feature(ltmp[2], ltmp[3])]

										ltmp_sets.append(["|".join(iset), dpaths[lpid[d[vid][n]]]])

								if len(ltmp_sets) > 0:
									lsets.append(ltmp_sets)

	return lsets



# --
# find_alt_starts
# 
def find_alt_starts(g, dtid):

	lsets = []
	all_introns = False
	dpaths = {}

	# get the chrom name
	ntmp = g.vs[0]['name']
	chrom = ntmp.split(":")[0]

	# get all output lists - first is the list, per node, of adjacent nodes
	# and the second is the same but the inner list elements are edges

	lvout = g.get_adjlist(mode=ig.OUT)
	lvin = g.get_adjlist(mode=ig.IN)
	leout = g.get_inclist(mode=ig.OUT)
	lein = g.get_inclist(mode=ig.IN)

	for i in range(len(lvin)):
		if len(lvin[i]) > 1:
			# this node has multiple inputs - are they all introns?
			all_introns = True
			for j in lein[i]:
				if g.es[j]['type'] != "intron":
					all_introns = False
					
			if not all_introns:
				continue

			# tid set
			ltid = g.vs[i]['tid_set']
			name0 = g.vs[i]['name']

			# build a dict keyed by the paths found
			dpaths = {}
			# per transcript id find the two step path along inputs by vertex names. 
			# confirm that 
			# the last of the two steps is the last step in the transcript
			for tid in ltid:
				sz_path = name0
				vid_set = dtid[tid]
				g_tid = g.subgraph(vid_set)
				# remove any edges that do not belong to this tid
				edrops = []
				for k in range(len(g_tid.es)):
					if tid not in g_tid.es[k]['tid_set']:
						edrops.append(k)

				# drop edges
				if len(edrops) > 0:
					g_tid.delete_edges(edrops)

				# find the first node
				lvidx = [g_tid.vs.find(name=name0).index]
				g_tid_lin = g_tid.get_adjlist(mode=ig.IN)
				
				for j in range(2):
					# get next and next
					if len(g_tid_lin[lvidx[j]]) > 0:
						lvidx.append(g_tid_lin[lvidx[j]][0])
					else:
						break

				# confirm two steps and the last one is the end of the 
				# transcript

				if len(lvidx)==3:
					if len(g_tid_lin[lvidx[2]])==0:
						pid = "|".join(g_tid.vs[lvidx]['name'])
						if pid not in dpaths:
							dpaths[pid] = []

						dpaths[pid].append(tid)



			lpid = dpaths.keys()
			if len(lpid) > 1:
				
				# there is more than one end from the initial node. 
				# we need a set of invalid pairs. those would be ones with
				# identical introns or matching end nodes

				# use a graph to figure this out by finding all of the pairs
				# that are valid and then clustering them. the clustering has 
				# to require that within each cluster all nodes are connected to 
				# all other nodes (i.e. all pairs of nodes within a cluster must be
				# a valid pair)

				# first break up each of the paths into node names
				llpid = []
				for j in range(len(lpid)):
					ltmp = lpid[j].split("|")
					llpid.append(ltmp)


				# start the graph
				gpaths = ig.Graph()
				gpaths.add_vertices(len(lpid))
				gedges = []

				# look for valid pairs
				for j in range(len(lpid)-1):
					for k in range(j+1, len(lpid)):
						# valid pair if they don't share second or third
						# positions
						if llpid[j][1] != llpid[k][1] and llpid[j][2] != llpid[k][2]:
							# valid pair
							gedges.append((j, k))

				if len(gedges) > 0:
					# add the edges to the graph
					gpaths.add_edges(gedges)
					# cluster
					groups = gpaths.clusters(mode=ig.STRONG).membership

					# evaluate the clusters
					lpid_sets = {}

					# based on the clustering, bundle the lpid indices into sets
					for j in range(len(groups)):
						if groups[j] not in lpid_sets:
							lpid_sets[groups[j]] = []

						lpid_sets[groups[j]].append(j)

					# if there are any sets now we can generate the results

					for j in lpid_sets.keys():
						ltmp_sets = []
						if len(lpid_sets[j]) > 1:
							# we have a set
							for k in range(len(lpid_sets[j])):
								# make the intron id
								ltmp = lpid[lpid_sets[j][k]].split("|")
								iid = combine_vid_to_feature(ltmp[0], ltmp[1])
								ltmp_sets.append([iid, dpaths[lpid[lpid_sets[j][k]]]])

							# scan these to see if any of them have identical introns so they can 
							# be collapsed
							d = {}
							for k in range(len(ltmp_sets)):
								if ltmp_sets[k][0] not in d:
									d[ltmp_sets[k][0]] = set([])

								d[ltmp_sets[k][0]].update(ltmp_sets[k][1])

							if len(d.keys()) < len(ltmp_sets):
								ltmp_sets = []
								for iid in d.keys():
									ltmp_sets.append([iid, list(d[iid])])

							lsets.append(ltmp_sets)
						
	return lsets
# --
# find_alt_ends
# 
def find_alt_ends(g, dtid):

	lsets = []
	all_introns = False
	dpaths = {}

	# get the chrom name
	ntmp = g.vs[0]['name']
	chrom = ntmp.split(":")[0]

	# get all output lists - first is the list, per node, of adjacent nodes
	# and the second is the same but the inner list elements are edges

	lvout = g.get_adjlist(mode=ig.OUT)
	lvin = g.get_adjlist(mode=ig.IN)
	leout = g.get_inclist(mode=ig.OUT)
	lein = g.get_inclist(mode=ig.IN)

	for i in range(len(lvout)):
		if len(lvout[i]) > 1:
			# this node has multiple outputs - are they all introns?
			all_introns = True
			for j in leout[i]:
				if g.es[j]['type'] != "intron":
					all_introns = False
					
			if not all_introns:
				continue

			# tid set
			ltid = g.vs[i]['tid_set']
			name0 = g.vs[i]['name']

			# build a dict keyed by the paths found
			dpaths = {}
			# per transcript id find the two step path by vertex names. confirm that 
			# the last of the two steps is the last step in the transcript
			for tid in ltid:
				sz_path = name0
				vid_set = dtid[tid]
				g_tid = g.subgraph(vid_set)
				# remove any edges that do not belong to this tid
				edrops = []
				for k in range(len(g_tid.es)):
					if tid not in g_tid.es[k]['tid_set']:
						edrops.append(k)

				# drop edges
				if len(edrops) > 0:
					g_tid.delete_edges(edrops)

				# find the first node
				lvidx = [g_tid.vs.find(name=name0).index]
				g_tid_lout = g_tid.get_adjlist(mode=ig.OUT)
				
				for j in range(2):
					# get next and next
					if len(g_tid_lout[lvidx[j]]) > 0:
						lvidx.append(g_tid_lout[lvidx[j]][0])
					else:
						break

				# confirm two steps and the last one is the end of the 
				# transcript

				if len(lvidx)==3:
					if len(g_tid_lout[lvidx[2]])==0:
						pid = "|".join(g_tid.vs[lvidx]['name'])
						if pid not in dpaths:
							dpaths[pid] = []

						dpaths[pid].append(tid)



			lpid = dpaths.keys()
			if len(lpid) > 1:
				
				# there is more than one end from the initial node. 
				# we need a set of invalid pairs. those would be ones with
				# identical introns or matching end nodes

				# use a graph to figure this out by finding all of the pairs
				# that are valid and then clustering them. the clustering has 
				# to require that within each cluster all nodes are connected to 
				# all other nodes (i.e. all pairs of nodes within a cluster must be
				# a valid pair)

				# first break up each of the paths into node names
				llpid = []
				for j in range(len(lpid)):
					ltmp = lpid[j].split("|")
					llpid.append(ltmp)

				# start the graph
				gpaths = ig.Graph()
				gpaths.add_vertices(len(lpid))
				gedges = []

				# look for valid pairs
				for j in range(len(lpid)-1):
					for k in range(j+1, len(lpid)):
						# valid pair if they don't share second or third
						# positions
						if llpid[j][1] != llpid[k][1] and llpid[j][2] != llpid[k][2]:
							# valid pair
							gedges.append((j, k))

				if len(gedges) > 0:
					# add the edges to the graph
					gpaths.add_edges(gedges)
					# cluster
					groups = gpaths.clusters(mode=ig.STRONG).membership


					# evaluate the clusters
					lpid_sets = {}

					# based on the clustering, bundle the lpid indices into sets
					for j in range(len(groups)):
						if groups[j] not in lpid_sets:
							lpid_sets[groups[j]] = []

						lpid_sets[groups[j]].append(j)

					# if there are any sets now we can generate the results

					for j in lpid_sets.keys():
						ltmp_sets = []
						if len(lpid_sets[j]) > 1:
							# we have a set
							for k in range(len(lpid_sets[j])):
								# make the intron id
								ltmp = lpid[lpid_sets[j][k]].split("|")
								iid = combine_vid_to_feature(ltmp[0], ltmp[1])
								ltmp_sets.append([iid, dpaths[lpid[lpid_sets[j][k]]]])

							# scan these to see if any of them have identical introns so they can 
							# be collapsed
							d = {}
							for k in range(len(ltmp_sets)):
								if ltmp_sets[k][0] not in d:
									d[ltmp_sets[k][0]] = set([])

								d[ltmp_sets[k][0]].update(ltmp_sets[k][1])

							if len(d.keys()) < len(ltmp_sets):
								ltmp_sets = []
								for iid in d.keys():
									ltmp_sets.append([iid, list(d[iid])])

							lsets.append(ltmp_sets)						
	return lsets


# --
# load_gtf
# loads the GTF row by row into a dict keyed by chromosome and then by 
# transcript id
def load_gtf(f):

	dchroms = {}
	szl = ""
	aln = []
	tattr = []
	has_tattr = False
	attr = []

	fin = open(f, "r")

	for szl in fin:
		aln = szl.strip().split("\t")

		if aln[_TYPE] == "transcript":
			tattr = parse_gtf_attr(aln[_ATTR])
			has_tattr = True
			continue

		if aln[_TYPE] == "exon":

			attr = parse_gtf_attr(aln[_ATTR])

			# add to dict
			if aln[_CHROM] not in dchroms:
				# add new dict for this chromosome
				dchroms[aln[_CHROM]] = {}

			if attr['transcript_id'] not in dchroms[aln[_CHROM]]:
				# add new list for this transcript
				dchroms[aln[_CHROM]][attr['transcript_id']] = []

			if has_tattr:
				aln.append(tattr)
			else:
				aln.append(attr)
				
			# add this line
			dchroms[aln[_CHROM]][attr['transcript_id']].append(list(aln))

	fin.close()

	return dchroms

# --
# parse_gtf_attr
# parses the attributes field of a GTF row into a dict keyed by the 
# available attribute names
def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	#
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

# --
# transcript_to_da
# extract the donor/acceptor sites as well as exon and introns from a 
# transcript
def transcript_to_ve(lt, dnodes, eset, iset):

	pstrand = get_strand(lt[0])=="+"
	lexons = []
	lintrons = []

	tid = get_tid(lt[0])

	# confirm the rows are sorted
	estarts = []
	for i in range(len(lt)):
		estarts.append(get_start(lt[i]))

	o = order(estarts)

	i = 0
	vidx = len(dnodes.keys())

	while i < len(o):
		# add nodes - start and end of feature

		eid = ""
		iid = ""

		vid = str(get_start(lt[o[i]]))
		if vid not in dnodes:
			dnodes[vid] = [vidx, []]
			vidx += 1

		dnodes[vid][1].append(tid)

		vid = str(get_end(lt[o[i]]))
		if vid not in dnodes:
			dnodes[vid] = [vidx, []]
			vidx += 1

		dnodes[vid][1].append(tid)

		# add exon edge
		if pstrand:
			eid = ":".join(map(str, [get_start(lt[o[i]]), get_end(lt[o[i]])]))
		else:
			eid = ":".join(map(str, [get_end(lt[o[i]]), get_start(lt[o[i]])]))

		if eid not in eset:
			eset[eid] = []

		eset[eid].append(tid)

		if i > 0:
			# add intron
			if pstrand:
				iid = ":".join(map(str, [get_end(lt[o[i-1]]), get_start(lt[o[i]])]))
			else:
				iid = ":".join(map(str, [get_start(lt[o[i]]), get_end(lt[o[i-1]])]))

			if iid not in iset:
				iset[iid] = []

			iset[iid].append(tid)

		# increment!
		i += 1


	# done, dnodes, eset and iset are modified in the calling code
	return 0



# --
# feature_overlap
# Returns True if the features overlap.  The features are specified in the format
# chrom:start-end
def feature_overlap(a, b):
	# explode a
	la = map(int, re.sub("^[^\:]+\:", "", a).split("-"))
	lb = map(int, re.sub("^[^\:]+\:", "", b).split("-"))

	if la[1] >= lb[0] and lb[1] >= la[0]:
		return True	

	return False

# --
# combine_vid_to_feature
# merge two node ids (format: chrom:position) into a feature id (chrom:start-end)
def combine_vid_to_feature(v1, v2):
	v1s = v1.split(":")
	v2s = v2.split(":")

	fid = "{}-{}".format(v1, v2s[1])

	if int(v2s[1]) < int(v1s[1]):
		# the positions are reversed
		fid = "{}-{}".format(v2, v1s[1])

	return fid

# --
# trim_feature
# this function is used to trim a single base off of both ends of the
# feature to make it compatible with the coordinates of an intron rather
# that the end positions of an exon
def trim_feature(v):
	s1 = v.split(":")
	s2 = s1[1].split("-")

	vHat = "{}:{}-{}".format(s1[0], int(s2[0])+1, int(s2[1])-1)
	return(vHat)


def order(v):
	return np.argsort(v)

# helpers

def get_chrom(lgrow):
	return(lgrow[0])

def get_chrom(lgrow):
	return(lgrow[2])

def get_start(lgrow):
	return(int(lgrow[3]))

def get_end(lgrow):
	return(int(lgrow[4]))

def get_strand(lgrow):
	return(lgrow[6])

def get_tid(lgrow):
	if len(lgrow) > 9:
		return(lgrow[9]['transcript_id'])

	attr = parse_gtf_attr(lgrow[8])
	return(attr['transcript_id'])

def get_gid(lgrow):
	if len(lgrow) > 9:
		return(lgrow[9]['gene_id'])

	attr = parse_gtf_attr(lgrow[8])
	return(attr['gene_id'])

def get_gname(lgrow):
	res = "NA"
	if len(lgrow) > 9:
		if "gene_name" in lgrow[9]:
			res = lgrow[9]['gene_name']

	attr = parse_gtf_attr(lgrow[8])
	if "gene_name" in attr:
		res = attr['gene_name']

	return res


# --
# build_chrom_graph
# builds graph from a single chromosome of transcripts
def build_chrom_graph(didx):

	return 0

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Build alternative splicing index.")
parser.add_argument('gtf', type=str, 
	help="GTF to base index on.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

