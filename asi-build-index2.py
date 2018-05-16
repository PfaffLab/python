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
# Same as asi-build-index.py but this one creates a GFF type annotation such 
# as those provided with MISO describing genes and transcripts that make up 
# alternative paths.  Output provides a "gene" row and then "mRNA" rows for 
# each overall path and finally the exon rows that belong to each path.
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser, basename

from igraph import *

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

	# load the GTF
	sys.stderr.write("Loading {}\n".format(args.gtf))
	dannot = parse_gtf(args.gtf)
	sys.stderr.write("sorting exons within isoforms\n")
	parsed_gtf_sort_exons(dannot)

	# by bundling 5' and 3' ends of junctions we can search out all of the 
	# alternative splice sites and figure out what's going on fairly quickly

	# build and exon, intron, 5p and 3p database. 

	dexons = {}
	dintrons = {}
	d5p = {}
	d3p = {}

	lverts = []
	ledges = []

	eid = ""
	eid_last = ""
	iid = ""
	p5id = ""
	p3id = ""

	# these dicts are used to link 5p to 3p and 3p to 5p

	sys.stderr.write("assembling 5p and 3p dicts\n")

	d5to3 = {}
	d3to5 = {}

	for tid in dannot.keys():

		for i in range(len(dannot[tid])):

			eid = "{}:{}:{}".format(dannot[tid][i].rname, 
				dannot[tid][i].start, dannot[tid][i].end)

#			if eid not in dexons:
#				dexons[eid] = []

			# append in the transcript, index and length
#			dexons[eid].append((tid, i, dannot[tid][i].end-dannot[tid][i].start+1))

			if i > 0:
				# make the intron
				iid = "{}:{}:{}".format(dannot[tid][i].rname, 
					dannot[tid][i-1].end+1, dannot[tid][i].start-1)

#				if iid not in dintrons:
#					dintrons[iid] = []

				# append the transcript and index but also the upstream and downstream exons
#				dintrons[iid].append((tid, i-1, eid_last, eid))

				# if we have an intron then we also have a 5p and 3p end

				p5id = "{}:{}:5p".format(dannot[tid][i].rname, dannot[tid][i-1].end+1)
				if p5id not in d5p:
					d5p[p5id] = []

				# append in the transcript and the upstream exon and downstream intron
				d5p[p5id].append((tid, eid_last, iid))

				p3id = "{}:{}:3p".format(dannot[tid][i].rname, dannot[tid][i].start-1)
				if p3id not in d3p:
					d3p[p3id] = []

				# append in the transcript and the upstream exon and downstream intron
				d3p[p3id].append((tid, eid, iid))

#				if p5id not in d5to3:
#					d5to3[p5id] = set()
#				
#				d5to3[p5id].update([p3id])
#
#				if p3id not in d3to5:
#					d3to5[p3id] = set()
#
#				d3to5[p3id].update([p5id])

			eid_last = eid

	# using the 5p and 3p dicts we can make a graph of transcript to transcript connections
	# based on sharing of donor/acceptor sites. then we can 

	sys.stderr.write("building graph\n")

	dverts = {}
	ledges = []
	lvname = []
	lvid = []
	vidx = 0

	for p5id in d5p.keys():
		if len(d5p[p5id]) > 1:
			# shared donor site. link all of the transcript ids as edges
			lvid = []
			for i in range(len(d5p[p5id])):
				tid = d5p[p5id][i][0]
				if tid not in dverts:
					dverts[tid] = vidx
					lvname.append(tid)
					vidx += 1

				lvid.append(dverts[tid])

			# make pairs of the vertex ids and add as edges
			for i in range(len(lvid)-1):
				for j in range(i+1, len(lvid)):
					ledges.append((lvid[i], lvid[j]))

	for p3id in d3p.keys():
		if len(d3p[p3id]) > 1:
			# shared acceptor site. link all of the transcript ids
			lvid = []
			for i in range(len(d3p[p3id])):
				tid = d3p[p3id][i][0]
				if tid not in dverts:
					dverts[tid] = vidx
					lvname.append(tid)
					vidx += 1

				lvid.append(dverts[tid])

			# make pairs of the vertex ids and add as edges
			for i in range(len(lvid)-1):
				for j in range(i+1, len(lvid)):
					ledges.append((lvid[i], lvid[j]))

	#
	# build graph to cluster the transcript ids
	#

	g = Graph(directed=False)
	g.add_vertices(len(dverts.keys()))
	g.add_edges(ledges)

	##
	# cluster the graph (i.e. cluster transcripts by common donor/acceptor connections)
	sys.stderr.write("clustering transcripts in graph\n")
	g_clusters = g.clusters(mode=WEAK)
	v_groups = g_clusters.membership

	##
	# get group indices
	tid_sets = {}
	for i in range(len(v_groups)):
		setid = "TSET_{:08d}".format(v_groups[i])

		if setid not in tid_sets:
			tid_sets[setid] = []

		tid_sets[setid].append(lvname[i])

	sys.stderr.write("searching for alternative splice paths in transcript clusters\n")

	lSE = [] # single exon skip (alt cassette)
	lTE = [] # multi exon (tandem cassette)
	lMX = [] # mutually exclusive
	lAS = [] # alt start
	lAE = [] # alt end
	lRI = [] # retained intron
	lA5 = [] # alt 5p
	lA3 = [] # alt 3p

	dalts = {} # dict to manage the alt structures. this will help bundle alt pairs

	#
	# now loop through the transcript clusters and parse out exons and introns	
	for tsetid in tid_sets:

		##
		# find the largest number of exons in any one transcript in this cluster
		max_exons = 0
		multi_iso = len(tid_sets[tsetid]) > 1

		for tid in tid_sets[tsetid]:
			if len(dannot[tid]) > max_exons:
				max_exons = len(dannot[tid])

		##
		#
		# cassette search
		#
		##

		if max_exons > 2 and multi_iso:

			# maximum number of skipped exons would be the longest number possible 
			# minus 2
			max_skips = max_exons-2

			#
			# build a list of each type of combination of exons (no skip, 1 skipped, etc)
			# and store them by path id which is the [=====A--------B=====] A:B coordinates
			# of the inner edges of the first and last exon of the path. 
			#

			# setup the list of dicts necessary for all possible exon skip combinations
			esets = []
			for i in range(max_skips+1):
				esets.append({})

			# populate the dicts. also build an exon dict to track what transcripts each 
			# exon belongs to

			dexons = {}
			for tid in tid_sets[tsetid]:
				num_exons = len(dannot[tid])
				
				# add the first exon into the exons dict
				eid = dannot[tid][0].feature_id()
				if eid not in dexons:
					dexons[eid] = set()

				dexons[eid].update([tid])

				for i in range(num_exons-1):
					for j in range(i+1, num_exons):

						# add this exon if necessary
						eid = dannot[tid][j].feature_id()
						if eid not in dexons:
							dexons[eid] = set()

						dexons[eid].update([tid])

						# get skip size
						skip_size = (j - i) - 1

						# make the path id
						path_id = "{}:{}:{}".format(dannot[tid][i].rname, dannot[tid][i].end, dannot[tid][j].start)

						if path_id not in esets[skip_size]:
							# at each path i want the transcript id and the indices of the exons
							esets[skip_size][path_id] = []

						esets[skip_size][path_id].append((tid, i, j))

			#==================================================================
			#
			# find alt cassette and tandem cassette features
			#
			#==================================================================			

			#
			# now we can use these dicts to search for cassette and mutually exclusive events.
			# 

			# check for single skips
			lpaths = list(set(esets[0].keys()).intersection(esets[1].keys()))
			if len(lpaths) > 0:
				# got some
				for path_id in lpaths:
					
					if path_id not in dalts:
						dalts[path_id] = []

					dalts[path_id].append(dict(type="SE", id=path_id, a=esets[0][path_id], b=esets[1][path_id]))


			# look for tandem skips
			for i in range(2, len(esets)):
				lpaths = list(set(esets[0].keys()).intersection(esets[i].keys()))
				if len(lpaths) > 0:
					for path_id in lpaths:
						if path_id not in dalts:
							dalts[path_id] = []

						dalts[path_id].append(dict(type= "TE", id=path_id, a=esets[0][path_id], b=esets[i][path_id]))

			#==================================================================
			#
			# find mutually exclusive features
			#
			#==================================================================			

			# look for mutually exclusive.  these will be path ids in esets[1] that have more than
			# one transcript. we will have to fetch the coordinates of the skipped exon and see 
			# if there are alternative cases.  we also have to know if those two exons ever appear
			# in the same isoform together, though

			for path_id in esets[1].keys():
				if len(esets[1][path_id]) > 1:
					# for each of these paths we need to fetch the middle exon
					lmexon = set()
					eid2p = {}
					i = 0
					for p in esets[1][path_id]:
						eid = dannot[p[0]][p[1]+1].feature_id()
						lmexon.update([eid])
						if eid not in eid2p:
							eid2p[eid] = []

						eid2p[eid].append(i)
						i += 1

					if len(lmexon) > 1:

						# i guess at this point each one of these has to be compared to 
						# the other to see if we have a mutually exclusive set. that means
						# they cannot share any transcript ids and they cannot share 
						# their 5p or 3p ends
						
						lmexon = list(lmexon)
						n = len(lmexon)
						
						for i in range(n-1):
							for j in range(i+1, n):
								# assume we found one then check for reasons that it 
								# would not be a mutually exclusive set
								bpass = True
								leid1 = lmexon[i].split(":")
								leid2 = lmexon[j].split(":")

								if leid1[1]==leid2[1] or leid1[2]==leid2[2]:
									bpass = False

								tid_int = set(dexons[lmexon[i]]).intersection(lmexon[j])
								if len(tid_int) > 0:
									bpass = False

								if bpass:
									# got it
									pa = []
									for k in eid2p[lmexon[i]]:
										pa.append(esets[1][path_id][k])
									pb = []
									for k in eid2p[lmexon[j]]:
										pb.append(esets[1][path_id][k])

									if path_id not in dalts:
										dalts[path_id] = []
									
									dalts[path_id].append(dict(type="MX", id=path_id, a=pa, b=pb))

		# -- end alt cassette and mutually exlusive search

		#==================================================================
		#
		# find retained intron and alt 3'/5' features
		#
		#==================================================================			

		# 
		# to find retained intron and alt 3' and 5' sets we have to build a new set
		# of combinations that use the following anchors A and B: A===]------[===B .
		# before we had exon skips as a factor but this time we'll do each exon and then 
		# each pair of exons.  for intron retention we'll check between exons and pairs
		# and for alt 3' and alt 5' we'll check within the pairs
		#

		# retained intron requires there to be at least two isoforms and at least one 
		# with more than one exon. same goes for alt 3' and alt 5'
		if max_exons > 1 and multi_iso:

			esets = [{}, {}]

			for tid in tid_sets[tsetid]:
				num_exons = len(dannot[tid])
				
				for i in range(num_exons):
					# put the exon in as a path
					path_id = dannot[tid][i].feature_id()
					if path_id not in esets[0]:
						esets[0][path_id] = []

					esets[0][path_id].append((tid, i, -1))

					if i > 1:
						# add to esets[1]

						path_id = "{}:{}:{}".format(dannot[tid][i].rname, 
							dannot[tid][i-1].start, dannot[tid][i].end)

						if path_id not in esets[1]:
							esets[1][path_id] = []

						esets[1][path_id].append((tid, i-1, i))

			# retained introns will be matches between esets[0] and esets[1]
			lpaths = set(esets[0].keys()).intersection(esets[1].keys())
			if len(lpaths) > 0:
				# got some
				for path_id in lpaths:
					if path_id not in dalts:
						dalts[path_id] = []

					dalts[path_id].append(dict(type="RI", id=path_id, a=esets[0][path_id], b=esets[1][path_id]))

			# to find alt 3' and alt 5' we have to check for paths in esets[1] that
			# have multiple entries and at those entries we have to figure out if there
			# are multiple paths
			for path_id in esets[1].keys():
				if len(esets[1][path_id]) > 1:
					#print esets[1][path_id]
					dpaths = {}
					for i in range(len(esets[1][path_id])):
						p = esets[1][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths:
							dpaths[pid] = []

						dpaths[pid].append(esets[1][path_id][i])

					# if there is only 1 path then we're done with this, otherwise
					# we do have something or possibly multiple things
					kpaths = dpaths.keys()
					if len(kpaths) > 1:
						for i in range(len(kpaths)-1):
							for j in range(i+1, len(kpaths)):
								eidA = kpaths[i].split("-")
								eidB = kpaths[j].split("-")

								if eidA[0] == eidB[0] and eidA[1] != eidB[1]:
									# alt 3'
									if path_id not in dalts:
										dalts[path_id] = []

									dalts[path_id].append(dict(type="A3", id=path_id, 
										a=dpaths[kpaths[i]], b=dpaths[kpaths[j]]))

								elif eidA[0] != eidB[0] and eidA[1] == eidB[1]:
									# alt 5'
									if path_id not in dalts:
										dalts[path_id] = []

									dalts[path_id].append(dict(type="A5", id=path_id, 
										a=dpaths[kpaths[i]], b=dpaths[kpaths[j]]))
		
		# -- end retained intron, alt 3' and 5' search

		#==================================================================
		#
		# find alt starts
		#
		#==================================================================			

		# collect all exon 1-2 pairs. id will be the acceptor of the second
		# exon. then in a second set collect all pairs that are NOT 1-2

		if multi_iso and max_exons > 1:

			esets = [{}, {}]

			for tid in tid_sets[tsetid]:
				num_exons = len(dannot[tid])
				
				for i in range(1, num_exons):
					if i > 1:
						# add to esets[1]
						path_id = "{}:{}".format(dannot[tid][i].rname, 
							dannot[tid][i].start)

						if path_id not in esets[1]:
							esets[1][path_id] = []

						esets[1][path_id].append((tid, i-1, i))

					else:
						# add to esets[0]
						path_id = "{}:{}".format(dannot[tid][i].rname, 
							dannot[tid][i].start)

						if path_id not in esets[0]:
							esets[0][path_id] = []

						esets[0][path_id].append((tid, i-1, i))


			# check within the first set. this is sort of like above where we looked
			# for alt 3' and 5'
			for path_id in esets[0].keys():
				if len(esets[0][path_id]) > 1:
					
					# make a lookup to track the path ids					
					dpaths = {}
					for i in range(len(esets[0][path_id])):
						p = esets[0][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths:
							dpaths[pid] = []

						dpaths[pid].append(esets[0][path_id][i])


					# if there is only 1 path then we're done with this, otherwise
					# we do have something or possibly multiple things
					kpaths = dpaths.keys()
					if len(kpaths) > 1:
						for i in range(len(kpaths)-1):
							for j in range(i+1, len(kpaths)):
								# split the pids into the pair of exons in each case.
								# then since we are only interested in the end exon we have
								# to extract that feature's coordinates to compare between 
								# the two cases
								eidA = kpaths[i].split("-")
								eidA_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidA[0])[0]
								eidB = kpaths[j].split("-")
								eidB_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidB[0])[0]

								if eidA_edges[0] != eidB_edges[0] and eidA_edges[1] != eidB_edges[1]:
									# got it
									if path_id not in dalts:
										dalts[path_id] = []

									dalts[path_id].append(dict(type="AS", id=path_id, 
										a=dpaths[kpaths[i]], b=dpaths[kpaths[j]]))

			# check for overlap between first set and second set. 
			lpaths = set(esets[0].keys()).intersection(esets[1].keys())

			if len(lpaths) > 0:
				for path_id in lpaths:
					# reference is the set from eset[0]. we need to check each unique 
					# exon pair in esets[0] against each in esets[1]
					dpaths_ref = {}
					dpaths_alt = {}

					for i in range(len(esets[0][path_id])):
						p = esets[0][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths_ref:
							dpaths_ref[pid] = []

						dpaths_ref[pid].append(esets[0][path_id][i])

					for i in range(len(esets[1][path_id])):
						p = esets[1][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths_alt:
							dpaths_alt[pid] = []

						dpaths_alt[pid].append(esets[1][path_id][i])

					# now check for pairs of first set to second set that do not
					# match
					for rpid in dpaths_ref.keys():
						for apid in dpaths_alt.keys():
							if rpid != apid:

								eidA = rpid.split("-")
								eidA_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidA[0])[0]
								eidB = apid.split("-")
								eidB_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidB[0])[0]

								if eidA_edges[0] != eidB_edges[0] and eidA_edges[1] != eidB_edges[1]:
									# got it
									if path_id not in dalts:
										dalts[path_id] = []

									dalts[path_id].append(dict(type="AS", id=path_id, 
										a=dpaths_ref[rpid], b=dpaths_alt[apid]))

			#==================================================================
			#
			# find alt ends
			#
			#==================================================================			

			# same idea as above but at the end of the transcripts instead

			esets = [{}, {}]

			for tid in tid_sets[tsetid]:
				num_exons = len(dannot[tid])
				
				for i in range(1, num_exons):
					if i < (num_exons-1):
						# add to esets[1]
						path_id = "{}:{}".format(dannot[tid][i].rname, 
							dannot[tid][i-1].end)

						if path_id not in esets[1]:
							esets[1][path_id] = []

						esets[1][path_id].append((tid, i-1, i))

					else:
						# add to esets[0]
						path_id = "{}:{}".format(dannot[tid][i].rname, 
							dannot[tid][i-1].end)

						if path_id not in esets[0]:
							esets[0][path_id] = []

						esets[0][path_id].append((tid, i-1, i))

			# search for combinations of ends that are both actual ends
			for path_id in esets[0].keys():
				if len(esets[0][path_id]) > 1:
					
					# make a lookup to track the path ids					
					dpaths = {}
					for i in range(len(esets[0][path_id])):
						p = esets[0][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths:
							dpaths[pid] = []

						dpaths[pid].append(esets[0][path_id][i])


					# if there is only 1 path then we're done with this, otherwise
					# we do have something or possibly multiple things
					kpaths = dpaths.keys()
					if len(kpaths) > 1:
						for i in range(len(kpaths)-1):
							for j in range(i+1, len(kpaths)):
								# split the pids into the pair of exons in each case.
								# then since we are only interested in the end exon we have
								# to extract that feature's coordinates to compare between 
								# the two cases
								eidA = kpaths[i].split("-")
								eidA_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidA[1])[0]
								eidB = kpaths[j].split("-")
								eidB_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidB[1])[0]

								if eidA_edges[0] != eidB_edges[0] and eidA_edges[1] != eidB_edges[1]:
									# got it
									if path_id not in dalts:
										dalts[path_id] = []

									dalts[path_id].append(dict(type="AE", id=path_id, 
										a=dpaths[kpaths[i]], b=dpaths[kpaths[j]]))

			# check for overlap between first set and second set. 
			lpaths = set(esets[0].keys()).intersection(esets[1].keys())

			if len(lpaths) > 0:
				for path_id in lpaths:
					# reference is the set from eset[0]. we need to check each unique 
					# exon pair in esets[0] against each in esets[1]
					dpaths_ref = {}
					dpaths_alt = {}

					for i in range(len(esets[0][path_id])):
						p = esets[0][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths_ref:
							dpaths_ref[pid] = []

						dpaths_ref[pid].append(esets[0][path_id][i])

					for i in range(len(esets[1][path_id])):
						p = esets[1][path_id][i]
						pid = "{}-{}".format(dannot[p[0]][p[1]].feature_id(), 
							dannot[p[0]][p[2]].feature_id())

						if pid not in dpaths_alt:
							dpaths_alt[pid] = []

						dpaths_alt[pid].append(esets[1][path_id][i])

					# now check for pairs of first set to second set that do not
					# match
					for rpid in dpaths_ref.keys():
						for apid in dpaths_alt.keys():
							if rpid != apid:

								eidA = rpid.split("-")
								eidA_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidA[1])[0]
								eidB = apid.split("-")
								eidB_edges = re.findall("^[^\:]+\:([0-9]+)\:([0-9]+)", eidB[1])[0]

								if eidA_edges[0] != eidB_edges[0] and eidA_edges[1] != eidB_edges[1]:
									# got it
									if path_id not in dalts:
										dalts[path_id] = []

									dalts[path_id].append(dict(type="AE", id=path_id, 
										a=dpaths_ref[rpid], b=dpaths_alt[apid]))


		# -- end alt start and alt end search

	# -- end main tid loop

	#
	# we need 2 things:
	#
	# 1. all of the paths that need to be quantified by alignments
	# 2. a table relating paths into pair-wise alt events
	#

	
	path_idx = 0
	path_id = ""
	# dict of paths based on distinct intron sets (or for retained introns just the exon)
	dpaths = {}
	# list of alt events by paired path ids (from dpaths) along with a type for the 
	# event
	devents = {}

	# going to put path to path connections into a graph
	dverts = {}
	lverts = []
	vidx = 0
	ledges = []
	ledge_types = []

	for pid in dalts.keys():
		print "# -- {}".format(pid)

		num_sets = len(dalts[pid])
		print num_sets

		# reassemble these in terms of intron sets
		for i in range(len(dalts[pid])):
			asets = dalts[pid][i]['a']
			bsets = dalts[pid][i]['b']

			daints = {}
			dbints = {}

			if dalts[pid][i]['type'] != "RI":

				for j in range(len(asets)):
					iid = get_full_path(dannot, asets[j][0], asets[j][1], asets[j][2])
					iid['set'] = asets[j]

					if iid['ipath'] not in daints:
						daints[iid['ipath']] = []

					daints[iid['ipath']].append(iid)

					if iid['ipath'] not in dpaths:
						dpaths[iid['ipath']] = []
						dverts[iid['ipath']] = vidx
						lverts.append(iid['ipath'])
						vidx += 1

					dpaths[iid['ipath']].append(iid)

				for j in range(len(bsets)):
					iid = get_full_path(dannot, bsets[j][0], bsets[j][1], bsets[j][2])
					iid['set'] = bsets[j]

					if iid['ipath'] not in dbints:
						dbints[iid['ipath']] = []

					dbints[iid['ipath']].append(iid)

					if iid['ipath'] not in dpaths:
						dpaths[iid['ipath']] = []
						dverts[iid['ipath']] = vidx
						lverts.append(iid['ipath'])
						vidx += 1

					dpaths[iid['ipath']].append(iid)

				for aid in daints.keys():
					for bid in dbints.keys():
						ledges.append((dverts[aid], dverts[bid]))
						ledge_types.append(dalts[pid][i]['type'])
	
	
	aidx = 0
	aid = ""
	aeid = ""
	aeidx = 0
	aid_paths = {} # dict of paths that should be quantified as a single gene id
	ae_sets = {}

	g = Graph(directed=False)
	g.add_vertices(len(lverts))
	g.add_edges(ledges)

	# annotate
	for i in range(len(lverts)):
		g.vs[i]['name'] = lverts[i]

	for i in range(len(ledges)):
		g.es[i]['type'] = ledge_types[i]
		g.es[i]['aid'] = lverts[ledges[i][0]]
		g.es[i]['bid'] = lverts[ledges[i][1]]

	g_clusters = g.clusters(mode=WEAK)
	v_groups = g_clusters.membership
	dgroups = split(range(len(lverts)), v_groups)

	for gid in dgroups.keys():
		ghat = g.subgraph(dgroups[gid])
		aid = "AID_{:08d}".format(aidx)
		aidx += 1

		aid_paths[aid] = []
		for i in range(len(ghat.vs)):
			aid_paths[aid].append(ghat.vs[i]['name'])
		


	return 0



#==============================================================================
# defs
#==============================================================================

def split(subject, factor):

	dout = {}
	# figure out unique factors
	for i in range(len(factor)):
		szf = "{}".format(factor[i])
		if szf not in dout:
			dout[szf] = []

		dout[szf].append(subject[i])

	return dout


#
# parse a GTF into a dict of GtfRow objects
def parse_gtf(fname):
	# variables
	gtfdb = {}
	grow = None

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:

		grow = GtfRow()
		grow.parse(szl.strip())

		if grow.type != "exon":
			continue

		if grow.tid not in gtfdb:
			gtfdb[grow.tid] = []

		gtfdb[grow.tid].append(grow)

	fin.close()

	return gtfdb

# 
# this should sort the passed dict in the calling code
def parsed_gtf_sort_exons(d):
	tid = ""
	for tid in d.keys():
		d[tid].sort(key=lambda x: x.start)

	return 0

#
# return a string with the full intron set coordinates from exon a to exon b
# in transcript tid
def get_intron_path(annot, tid, a, b):

	sz = ""

	if b > a:

		# start it
		sz = "{}:{}".format(annot[tid][a].rname, annot[tid][a].end)

		if b-a > 1:
			for i in range(a+1, b):
				sz = sz + "-{}:{}".format(annot[tid][i].start, annot[tid][i].end)

		# append last one
		sz = sz + "-{}".format(annot[tid][b].start)

	return sz

# 
# return a tuple with exon path string and list of exon lengths
def get_exon_path(annot, tid, a, b):

	elens = []
	
	sz = "{}:{}-{}".format(annot[tid][a].rname, annot[tid][a].start, annot[tid][a].end)
	elens.append(annot[tid][a].end-annot[tid][a].start+1)

	if b > a:
		for i in range(a+1, b+1):
			sz += ":{}-{}".format(annot[tid][i].start, annot[tid][i].end)
			elens.append(annot[tid][i].end-annot[tid][i].start+1)

	return (sz, elens)

def get_full_path(annot, tid, a, b):

	sz_intron = get_intron_path(annot, tid, a, b)
	sz_exon, elens = get_exon_path(annot, tid, a, b)

	chrom = annot[tid][0].rname
	bound = [annot[tid][a].start, annot[tid][b].end]

	return dict(ipath=sz_intron, epath=sz_exon, elens=elens, chrom=chrom, bound=bound)




#==============================================================================
# classes
#==============================================================================

#
# class to hold GTF row objects. 
class GtfRow(object):
	def __init__(self, rname="", db="", type="", start=0, end=0, strand="."):
		self.rname = rname
		self.db = db
		self.type = type
		self.start = start
		self.end = end
		self.strand = strand
		
		self.tid = ""
		self.gid = ""
		
		self.attrs = {}
		
	# parse info in GTF row, sz, into this object
	def parse(self, sz):
		aln = sz.strip().split("\t")
		self.rname = aln[0]
		self.db = aln[1]
		self.type = aln[2]
		self.start = int(aln[3])
		self.end = int(aln[4])
		self.strand = aln[6]
		
		# parse attributes from field 8
		fsplit = aln[8].split("\"")
		n = len(fsplit)-1
		i = 0
		while i < n:
			key = re.sub(';','',fsplit[i])
			self.attrs[key.strip()] = fsplit[i+1].strip()
			i += 2

		self.tid = self.attrs['transcript_id']
		self.gid = self.attrs['gene_id']
		
		return 0

	#--
	# tolist
	# return a list version of this with elements in place of the columns 
	# of a GTF
	def tolist(self):
		ltmp = [self.rname, self.db, self.type, 
			int(self.start), int(self.end), ".", self.strand, "."]
		
		szattr = "transcript_id \"{}\"; gene_id \"{}\";".format(self.tid, self.gid)
		akey = self.attrs.keys()
		if len(akey) > 0:
			# append additional attributes to the szattr string
			for aid in akey:
				szattr += " {} \"{}\";".format(aid, self.attrs[aid])
		
		ltmp.append(szattr)
		
		return(ltmp)

	def feature_id(self):
		sz = "{}:{}:{}".format(self.rname, self.start, self.end)
		return sz

	def feature_5p_id(self):
		sz = "{}:{}".format(self.rname, self.start)

	def feature_3p_id(self):
		sz = "{}:{}".format(self.rname, self.end)


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Build alternative splicing index GFF (like MISO style)")
parser.add_argument('gtf', type=str, 
	help="GTF to base index on.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

