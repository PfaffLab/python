#!/usr/bin/env python
#==============================================================================
# asi-build-indexG.py
#
# Shawn Driscoll
# 20160906
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Use graphs to find alternative splice patterns in genes.
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
# control the maximum number of skipped exons that may be returned in a
# tandem cassette event. this limit is in place mostly to improve speed
# with minor loss of splicing event information
MAX_SKIPS = 4

#==============================================================================
# main
#==============================================================================

def main(args):

	#--------------------------------------------------------------------------
	# 
	# define variables
	#
	#--------------------------------------------------------------------------

	#
	# use lists and dicts to assemble all of the information 
	# prior to actually adding nodes and edges to the graph. 
	# this goes much faster than building the graph as we go.
	tnodes = []
	tedges = []
	dtnodes = {}
	dedges = {}
	dverts = {}
	lverts = []
	ledges = []
	idx = 0
	vidx = 0
	eidx = 0
	tidx = 0
	num_verts = 0
	num_edges = 0
	num_paths = 0
	num_events = 0

	#
	# keep some running statistics on numbers of edges later when
	# building graphs for gene loci bundles
	ecount_mu = 0
	ecount_sigma = 0
	ecount0 = 0
	esigma0 = 0
	graph_count = 0
	num_tid = 0
	tid_set = None
	tid = None
	tidA = None
	tidB = None
	eid1 = None
	eid2 = None
	tid_groups = None
	tset_gid = None
	tset_gname = None


	# dict for all paths
	d_paths = {}
	d_literal_paths = {}
	# dict for all paired events by path id
	d_events = {}
	d_literal_events = {}
	# dict for all event loci
	d_locs = {}

	ae_id = ""
	ae_idx = 0
	aeloc_id = ""
	aeloc_idx = 0
	pe_id = ""
	pe_idx = 0

	tid_set = None
	in_nodes = None
	cid = None

	ltmp = None


	#--------------------------------------------------------------------------
	# 
	# get started
	#
	#--------------------------------------------------------------------------

	# load the GTF
	sys.stderr.write("Loading {}\n".format(args.gtf))
	dannot = parse_gtf(args.gtf)
	#
	# sort the exons within each transcript so that later when we 
	# build the graphs the exons are ordered and we won't have to 
	# check.
	sys.stderr.write("sorting exons within isoforms\n")
	parsed_gtf_sort_exons(dannot)


	# 
	# 20160912 - building graphs of entire chromosomes is taking too long so 
	# I'm going to bundle transcript ids up by shared 5' and 3' ends of exons
	# and then loop through all of those transcript bundles to build graphs.
	# that's also nice because if a bundle only has a single transcript then 
	# we're done with it and we can move on. 
	#

	sys.stderr.write("bundling transcript ids by shared donor/acceptor ends of exons\n")

	#
	# the following loop builds a single dict with 3'/5' id strings as keys
	# and at each index we will have a set of transcript ids that use the 
	# 3' or 5' site
	#

	dprimes = {}

	# loop through transcripts
	for tid in dannot.keys():
		# loop through exons in the transcript
		for exon in dannot[tid]:
			#
			# make 3' and 5' ids but append a 3/5 p identifyer because
			# i think there are some rare spots where a single position
			# is both due to anti-sense transcripts
			id3 = exon.feature_3p_id() + ":3p"
			id5 = exon.feature_5p_id() + ":5p"
			#
			# add in the 3' if not already there
			if id3 not in dprimes:
				# initalize empty set for transcript ids
				dprimes[id3] = set()
			#
			# add transcript id
			dprimes[id3].update([tid])

			#
			# add in the 5' if not already there
			if id5 not in dprimes:
				# initizlie empty set for transcript ids
				dprimes[id5] = set()
			#
			# add transcript id
			dprimes[id5].update([tid])

	#
	# finished collecing transcript ids by donor/acceptor ids. now build
	# a graph by adding transcript id nodes and edges connecting transcripts
	# that share any donor/acceptor sites. this basically will bundle transcripts
	# into loci
	#

	#
	# loop through 3/5' id keys
	for idprime in dprimes.keys():
		# is this necessary?
		#idx += 1
		#
		# count number of transcripts in this d/a list
		num_tid = len(dprimes[idprime])

		if num_tid < 2:
			# less than 2 means we can't make an edge so we don't need 
			# to add anything. this is also why later on the number of 
			# bundles is only equal to the number of bundles with alternative
			# splice events and not the total number of gene loci
			continue

		#
		# get list of transcript ids
		tid_set = list(dprimes[idprime])
		#
		# check if the tids are in already
		for tid in tid_set:
			if tid not in dtnodes:
				# add transcript id to the nodes list and the nodes dict
				tnodes.append(tid)
				# dict will hold the index of the transcript id that it 
				# will have later when added to the graph
				dtnodes[tid] = len(tnodes)-1

		# 
		# double loop through transcript id set to get all pairs of transcripts
		for i in range(num_tid-1):
			tidA = tid_set[i]
			for j in range(i+1, num_tid):
				tidB = tid_set[j]
				#
				# build ids for both possible orders of the transcript ids so we 
				# can check to see if we've already added this edge. this greatly
				# improves speed
				eid1 = "{}:{}".format(tidA, tidB)
				eid2 = "{}:{}".format(tidB, tidA)
				if (eid1 not in dedges) and (eid2 not in dedges):
					#
					# add the edge
					tedges.append((dtnodes[tidA], dtnodes[tidB]))
					# insert edge names into the edge dict so we don't add it again
					dedges[eid1] = 0
					dedges[eid2] = 0

	#
	# start graph
	tgraph = Graph(directed=False)
	#
	# add vertices - igraph only uses integers for vertices so we just 
	# add an anonymous number of nodes
	tgraph.add_vertices(len(tnodes))
	#
	# add edges
	tgraph.add_edges(tedges)
	#
	# cluster with weak connection (i.e. only a single link. i think this is only 
	# relevent when we have a directed graph anyways)
	tgraph_clust = tgraph.clusters(mode=WEAK)
	#
	# split up the transcript ids by the clustering memberships
	tid_groups = split(tnodes, tgraph_clust.membership)

	# let user know what we found
	sys.stderr.write("--> found {} transcript bundles\n".format(len(tid_groups)))

	#
	# save the locus bundles to a file
	sys.stderr.write("[saving bundle info to file]\n")
	fout = open("tbundles.log", "w")
	for tsetid in tid_groups.keys():
		# get transcript ids
		tid_set = tid_groups[tsetid]
		tset_gid = set()
		tset_gname = set()
		#
		# get gene names and ids from the parsed GTF
		for tid in tid_set:
			tset_gid.update([dannot[tid][0].attrs['gene_id']])
			tset_gname.update([dannot[tid][0].attrs['gene_name']])
		#
		# write a line for this transcript set
		fout.write(tsetid + "\t" + 
			",".join(list(tset_gid)) + "\t" + 
			",".join(list(tset_gname)) + "\t" + 
			",".join(tid_groups[tsetid]) + "\n")

	fout.close()

	# --
	# 
	# with all of the transcript ids clustered we can loop through those cluster sets
	# and build full graphs at each set then scan for alternative events
	#
	# --
	
	sys.stderr.write("scanning transcript bundles for alternative events\n")
	
	#
	# loop through transcript bundles
	for tsetid in tid_groups.keys():
		#
		# count and print progress message
		tidx += 1
		if (tidx % 1000) == 0:
			sys.stderr.write("processed {} bundles; {} paths in {} events in {} loci\n".format(tidx, num_paths, num_events, aeloc_idx))

		# get list of transcript ids
		tid_set = tid_groups[tsetid]

		# if less than 2 then we can skip it. interestingly this probably never
		# happens because of how we bundled the transcripts in the first place.
		if len(tid_set) < 2:
			continue

		# get chromosome id because I think I use it somewhere in this loop
		# from before when this iterated over chromosomes
		cid = dannot[tid_set[0]][0].rname

		#
		# graph these transcripts by building a graph of donor/acceptor nodes
		# and edges that are either intron or exon
		#

		dverts = {}
		lverts = []
		dedges = {}
		ledges = []
		vidx = 0
		eidx = 0
		num_verts = 0
		num_edges = 0

		if True:
			#
			# as before we are not going to build the graph as we go because
			# it slows execution. instead we'll build a dict of nodes and a
			# dict of edges then extract lists out of those to build the 
			# graph after this loop exits
			for tid in tid_set:

				# get exon count
				nexons = len(dannot[tid])

				# loop through exons
				for i in range(nexons):
					#
					# make acceptor/donor/exon id strings for use below	
					p5id, p3id, eid = exon_id_set(dannot[tid][i])

					#
					# check if acceptor id is in the dict
					if p5id not in dverts:
						# initalize acceptor with slots for transcript ids, name, type 
						# and index (future graph node index)
						dverts[p5id] = dict(name=p5id, tid=set(), type="a", index=vidx)
						# increment node counter
						vidx += 1

					#
					# check if donor id is in the dict
					if p3id not in dverts:
						# initalize as above with the acceptor id
						dverts[p3id] = dict(name=p3id, tid=set(), type="d", index=vidx)
						vidx += 1

					#
					# update transcript id sets for donor and acceptor nodes
					dverts[p5id]['tid'].update([tid])
					dverts[p3id]['tid'].update([tid])

					# 
					# check if this edge has been added already
					if eid not in dedges:
						# initalize edge with name, upstream acceptor and downstream donor, 
						# type (exons), transcript id set and the actual edge tuple that will 
						# be added to the graph later
						dedges[eid] = dict(name=eid, up=p5id, dn=p3id, type="e", exon5=0, exon3=0, 
							tid=set(), edge=(dverts[p5id]['index'], dverts[p3id]['index']), index=eidx)
						# increment edge counter
						eidx += 1

					#
					# update transcript id set for this edge
					dedges[eid]['tid'].update([tid])

					#
					# some exons may be used as starts for some isoforms but they'll be spliced
					# in in others so we need to flag those that are used as starts.
					# NOTE: we didn't do anything with this yet in terms of interpreting
					# alternative starts/ends. these would mark the type of exons that are
					# sometimes ends but sometimes fully spliced into a transcript.
					if i==0:
						dedges[eid]['exon5'] = 1
					elif i==(nexons-1):
						dedges[eid]['exon3'] = 1

					#
					# if we are at the second exon, at least, then we can add
					# the intron
					if i > 0:
						#
						# get donor/acceptor/intron id strings
						p5id, p3id, eid = intron_id_set(dannot[tid][i-1], dannot[tid][i])

						# 
						# check and initalize donor and acceptor ids if necessary. see above
						# for comments.
						if p5id not in dverts:
							dverts[p5id] = dict(name=p5id, tid=set(), type="d", index=vidx)
							vidx += 1

						if p3id not in dverts:
							dverts[p3id] = dict(name=p3id, tid=set(), type="a", index=vidx)
							vidx += 1

						#
						# update transcript id sets
						dverts[p5id]['tid'].update([tid])
						dverts[p3id]['tid'].update([tid])

						#
						# check if this edge has been added already
						if eid not in dedges:
							# initalize intron edge. see exon section for comments.
							dedges[eid] = dict(name=eid, up=p5id, dn=p3id, type="i", exon5=0, exon3=0, 
								tid=set(), edge=(dverts[p5id]['index'], dverts[p3id]['index']), index=eidx)
							eidx += 1

						#
						# update transcript set for this intron edge
						dedges[eid]['tid'].update([tid])

			# -- end tid loop

			#
			# finished assembling nodes and edges in dicts, now extract them into lists
			# and add them into a directed graph
			ghat = Graph(directed=True)
			#
			# add nodes
			vert_ids = dverts.keys()
			ghat.add_vertices(len(vert_ids))
			#
			# annotate nodes
			for vid in vert_ids:
				idx = dverts[vid]['index']
				ghat.vs[idx]['name'] = dverts[vid]['name']
				ghat.vs[idx]['tid'] = dverts[vid]['tid']
				ghat.vs[idx]['type'] = dverts[vid]['type']

			#
			# extract a list of the edge tuples. note by extracting the 
			# dict keys here we can make sure the edge annotation loop
			# executes in the same order as the edges are added.
			edge_ids = dedges.keys()
			ledges = []
			for eid in edge_ids:
				ledges.append(dedges[eid]['edge'])
			#
			# add all edges at once
			ghat.add_edges(ledges)
			#
			# loop through edge ids and annotate the edges
			eidx = 0
			for eid in edge_ids:
				ghat.es[eidx]['name'] = dedges[eid]['name']
				ghat.es[eidx]['up'] = dedges[eid]['up']
				ghat.es[eidx]['dn'] = dedges[eid]['dn']
				ghat.es[eidx]['type'] = dedges[eid]['type']
				ghat.es[eidx]['tid'] = dedges[eid]['tid']
				ghat.es[eidx]['exon5'] = dedges[eid]['exon5']
				ghat.es[eidx]['exon3'] = dedges[eid]['exon3']
				eidx += 1

		#
		# update edge count running stats
		tmp_count = ghat.ecount()
		graph_count += 1
		ecount_mu = ecount0 + (tmp_count-ecount0)/graph_count
		ecount_sigma = esigma0 + (tmp_count-ecount0)*(tmp_count-ecount_mu)
		ecount0 = ecount_mu
		esigma0 = ecount_sigma
		ecount_z = 0
		if graph_count > 2:
			ecount_z = (tmp_count-ecount0)/math.sqrt(esigma0/(graph_count-1))
		#
		# print a message if there is a graph with a crazy number of edges
		if ecount_z > 10:
			tmp_gname = set()
			for tid in tid_set:
				tmp_gname.update([dannot[tid][0].attrs['gene_name']])
			print tmp_gname
			sys.stderr.write("high edge count ({}): {}, {}\n".format(tsetid, tmp_count, ecount_z))

		# 
		# now scan graph for alternative events
		#


#			print "start exons"
#			e_starts = ghat.es.select(exon5_eq=1)
#			e_ends = ghat.es.select(exon3_eq=1)
#			for e in e_starts:
#				v = ghat.vs.select(name_eq=e['up'])
#				print e['name'], v['name'], v[0].index

		# make lists of alt donors, mode=out and mode=in
		# make lists of alt acceptors, mode=out and mode = in

		# 
		# get all donor nodes and all acceptor nodes
		vd = ghat.vs.select(type_eq="d")
		va = ghat.vs.select(type_eq="a")

		# 
		# initalize empty lists for the donors and acceptors we find 
		# to have more than 1 incoming or outgoing edge
		vd_aout = []
		vd_ain = []
		v_ends = []
		va_aout = []
		va_ain = []
		v_starts = []

		#
		# for each type of search (alt-donor to alt-acceptor, for example) we'll 
		# added in found paths to this list
		gid_path_sets = []

		gid_paths_found = {}

		#
		# loop through donors and acceptors, record node indices for those that 
		# have more than one incoming or outgoing edge.
		for v in vd:
			if v.degree(mode="OUT") > 1:
				vd_aout.append(v.index)
			elif v.degree(mode="OUT") == 0:
				# if node has no outgoing edge then it is an end node
				v_ends.append(v.index)

			if v.degree(mode="IN") > 1:
				vd_ain.append(v.index)

		for v in va:
			if v.degree(mode="OUT") > 1:
				va_aout.append(v.index)

			if v.degree(mode="IN") > 1:
				va_ain.append(v.index)
			elif v.degree(mode="IN") == 0:
				# if node has no incoming edge then it is a start
				v_starts.append(v.index)

		#
		# look for paths between donors with alt-out and acceptors with alt-in. this 
		# will reveal cassette and mutually exclusive types
		digout = dig_for_paths(ghat, vd_aout, va_ain)
		if len(digout) > 0:
			for dd in digout:
				gid_path_sets.append(dd)

		#
		# look for paths between acceptors with alt-out and acceptors with alt-in. this 
		# will reveal alt-5' types
		digout = dig_for_paths(ghat, va_aout, va_ain)
#		if digout is not None:
#			gid_path_sets.append(digout)
		if len(digout) > 0:
			for dd in digout:
				gid_path_sets.append(dd)


		#
		# look for paths between donors with alt-out and donors with alt-in. this will 
		# reveal alt-3' types
		digout = dig_for_paths(ghat, vd_aout, vd_ain)
#		if digout is not None:
#			gid_path_sets.append(digout)
		if len(digout) > 0:
			for dd in digout:
				gid_path_sets.append(dd)


		#
		# look for paths between acceptors with alt-out and donors with alt-in. this will 
		# reveal retained intron types.
		digout = dig_for_paths(ghat, va_aout, vd_ain)
#		if digout is not None:
#			gid_path_sets.append(digout)
		if len(digout) > 0:
			for dd in digout:
				gid_path_sets.append(dd)


		#
		# search for AS. paths from starts to alternative in's at acceptors that are 3 nodes long.
		# first we get all start to alt-in acceptor paths then narrow them down to 3 node length.
		# next we get the common alt-in acceptor anchors and do another path search for all 
		# 3 node paths incoming to those anchors. the returned paths are grouped by common 
		# introns and then paired off making sure that each pair contains at least 
		# one actual start exon. 
		digout = dig_for_paths(ghat, v_starts, va_ain, 1)		
		foo = pair_AS(ghat, digout)
		if len(foo) > 0:
			for dd in foo:
				gid_path_sets.append(dd)

		#
		# search for AE. the opposite search as that for AS.
		digout = dig_for_paths(ghat, vd_aout, v_ends, 1)		
		foo = pair_AE(ghat, digout)
		if len(foo) > 0:
			for dd in foo:
				gid_path_sets.append(dd)


		if False:

			#
			# look for alt start paths by exploring paths from alt start A ends of exons
			# to acceptor nodes with alt-ins. if any of the returned paths are 3 nodes long
			# then we can use the final edge of the path to explore paths going back upstream
			# from that edge
	
			aidx_used = []
			for didx in v_starts:
				for aidx in va_ain:
	
					if aidx in aidx_used:
						continue
	
					# make path transcript id set (common tids between start and end node)
					path_tid = ghat.vs[didx]['tid'].intersection(ghat.vs[aidx]['tid'])
					if len(path_tid) > 1:
						# this means more than one transcript id is shared which suggests there 
						# may be more than one intron/exon path between the nodes
						lpaths = find_all_paths(ghat, didx, aidx, path_tid)
						
						# make a set of valid start paths. these will be 3 nodes long and we'll 
						# use the end node of those sets to build up more sets
						valid_paths = []
						valid_path_tid = []
						for path in lpaths:
							if len(path) == 3:
								# this is actually the only stop we have to make because this reveals
								# that the current acceptor alt-in is an edge for an alt start event and 
								# since this back-tracks out of the edge to find all upstream paths we 
								# don't need to repeat it ever.
								
								if aidx in aidx_used:
									continue
	
								aidx_used.append(aidx)
								
								final_path_tid, path_valid = validate_path(ghat, path)
								
								if path_valid:
									
									# now we have to look at the end node and gather all of the 
									# valid two step paths back upstream
									
									back_paths0 = find_all_nstep_paths(ghat, path[2], 2, mode="IN")
									for p in back_paths0:
										p.reverse()
										btid, bvalid = validate_path(ghat, p)
										if bvalid:
											valid_paths.append(list(p))
											valid_path_tid.append(btid)
	
	
									#
									# one more step here. we have to make it so we only have a single
									# instance of each intron. the above code will pick up two paths for a 
									# single intron if the start exons are different lengths. 
									# example: 
									#   [==]-----------[====]
									# [====]-----------[====]
									# will come back as two paths. I want to keep the one of those with 
									# the longest start exon
									#
									if len(valid_paths) > 1:
										# make a new table
										dvalid_paths = {}
										for lp in valid_paths:
											#
											# setup an id based on the node indices of the 
											# donor and acceptor pair making up the intron
											did = "{}:{}".format(lp[1], lp[2])
											# 
											# get the length of the start exon by splicing it out of the
											# acceptor and donor names for that exon
											tmpA = ghat.vs[lp[0]]['name'].split(":")
											tmpB = ghat.vs[lp[1]]['name'].split(":")
											elen = float(tmpB[0])-float(tmpA[0])+1
											#
											# is it in the dict already?
											if did not in dvalid_paths:
												dvalid_paths[did] = dict(path=lp, elen=elen)
											else:
												# see if this elen is longer than the one we have already
												if elen > dvalid_paths[did]['elen']:
													dvalid_paths[did]['path'] = lp
													dvalid_paths[did]['elen'] = elen
	
										valid_paths = [dvalid_paths[did]['path'] for did in dvalid_paths.keys()]
	
									#
									# another thing we have check is that if the second node has alternative 
									# inputs (the donor end of the start exon) then they must all come from 
									# start nodes
									#
									if len(valid_paths) > 1:
										vtmp = []
										for lp in valid_paths:
											# check if this path is a start
											
											keep = True
	
											if lp[0] in v_starts:
												in_nodes = ghat.neighbors(lp[1], mode="IN")
												# if these are not all start nodes then we have to toss
												# this start because at least part of it is part of a spliced
												# exon which means this path's intron may end up being a part 
												# of some other feature like an SE event	
												for v in in_nodes:
													if v not in v_starts:
														keep = False
	
											if keep:
												vtmp.append(lp)
	
										valid_paths = list(vtmp)
	
									keep = False
									if len(valid_paths) > 1:
										# at least 2 paths
										# now we need to see if there is still a path that starts with a 
										# "start"
										for lp in valid_paths:
											if lp[0] in v_starts:
												keep = True
	
									# do we have more than one valid path?
									if (len(valid_paths) > 1) and keep:
										gid_path_sets.append(dict(p=valid_paths, t=valid_path_tid))
	
	
			#
			# this is the opposite idea from alt starts.  we have to find 3 node long
			# paths from donor alt outs to the ends
			#gid_path_sets = []
			didx_used = []
			for didx in vd_aout:
				if didx in didx_used:
					continue
	
				for aidx in v_ends:
					# make path transcript id set (common tids between start and end node)
					path_tid = ghat.vs[didx]['tid'].intersection(ghat.vs[aidx]['tid'])
					if len(path_tid) > 1:
						# this means more than one transcript id is shared which suggests there 
						# may be more than one intron/exon path between the nodes
						lpaths = find_all_paths(ghat, didx, aidx, path_tid)
						
						# make a set of valid start paths. these will be 3 nodes long and we'll 
						# use the end node of those sets to build up more sets
						valid_paths = []
						valid_path_tid = []
						for path in lpaths:
							if len(path) == 3:
	
								if didx in didx_used:
									continue
	
								didx_used.append(didx)
	
								final_path_tid, path_valid = validate_path(ghat, path)
								
								if path_valid:
									
									# since this is a valid 3node path from an alt donor to 
									# an actual end node now we want all 2 step (3 node) paths
									# from this same alt donor site..
									
									back_paths0 = find_all_nstep_paths(ghat, path[0], 2, mode="OUT")
									for p in back_paths0:
										# validate and append if it passes
										btid, bvalid = validate_path(ghat, p)
										if bvalid:
											valid_paths.append(list(p))
											valid_path_tid.append(btid)
	
									# rev
	#									valid_paths.append(path)
	#									valid_path_tid.append(final_path_tid)
	
									#
									# one more step here. we have to make it so we only have a single
									# instance of each intron. the above code will pick up two paths for a 
									# single intron if the end exons are different lengths. 
									# example: 
									# [====]-----------[==]
									# [====]-----------[=========]
									# will come back as two paths. I want to keep the one of those with 
									# the longest end exon
									#
									if len(valid_paths) > 1:
										# make a new table
										dvalid_paths = {}
										for lp in valid_paths:
											did = "{}:{}".format(lp[0], lp[1])
											tmpA = ghat.vs[lp[1]]['name'].split(":")
											tmpB = ghat.vs[lp[2]]['name'].split(":")
											elen = float(tmpB[0])-float(tmpA[0])+1
											#
											# is it in the dict already?
											if did not in dvalid_paths:
												dvalid_paths[did] = dict(path=lp, elen=elen)
											else:
												# see if this elen is longer than the one we have already
												if elen > dvalid_paths[did]['elen']:
													dvalid_paths[did]['path'] = lp
													dvalid_paths[did]['elen'] = elen
	
										valid_paths = [dvalid_paths[did]['path'] for did in dvalid_paths.keys()]
	
	
									#
									# opposite of the search we do above for starts.
									if len(valid_paths) > 1:
										vtmp = []
										for lp in valid_paths:
											#
											# check if this path is an end
											keep = True
											if lp[-1] in v_ends:
												out_nodes = ghat.neighbors(lp[1], mode="OUT")
												# keep it if all downstream edges are ends otherwise we have 
												# to toss it
												for v in out_nodes:
													if v not in v_starts:
														keep = False
	
											if keep:
												# keep it. if the path isn't a start then it's kept
												# by default
												vtmp.append(lp)
	
										valid_paths = list(vtmp)
	
									keep = False
									if len(valid_paths) > 1:
										# at least 2 paths
										# now we need to see if there is still a path that's actually an 'end'
										for lp in valid_paths:
											if lp[-1] in v_ends:
												keep = True
	
	
									# do we have more than one valid path?
									if (len(valid_paths) > 1) and keep:
										gid_path_sets.append(dict(p=valid_paths, t=valid_path_tid))

		#
		# if we don't have anything in gid_path_sets then we can move on to the next
		# transcript set
		if len(gid_path_sets) < 1:
			continue

		#
		# finished collecting bundles of valid paths, now the paths within each 
		# bundle are compared to discover alternative splice events.
		#

		# 
		# loop through bundles of paths
		for dbundle in gid_path_sets:
			# 
			# get paths found in the bundle
			valid_paths = dbundle['p']
			#
			# initialize some lists to contain the results
			bundle_results = []
			loc_ae_set = []
			loc_path_set = []

			#
			# evaluate all pairs of paths and determine if they form a
			# alternative event
			for i in range(len(valid_paths)-1):
				for j in range(i+1, len(valid_paths)):
					# compare paths
					rres, ref = compare_paths(ghat, valid_paths[i], valid_paths[j])
					# check result
					if rres is not None:
						# result is good so add the indices of the paths within
						# the 'valid_paths' list as well as the reference path
						# and the result (string indicating AE, AS, SE, TE..etc)
						bundle_results.append((i, j, ref, rres))

			#
			# check for results
			if len(bundle_results) > 0:
				#
				# something in the bundle made up an alternative event

				#
				# get all of the paths that are part of any event
				final_pid = set()					
				for i in range(len(bundle_results)):
					final_pid.update(bundle_results[i][0:2])

				# create path ids for each of the paths used
				path_ids = ["" for i in range(len(valid_paths))]
				# path "literal" ids which are build from the actual exon coodrinates in the path
				path_lids = list(path_ids)
				for i in list(final_pid):

					path_lids[i] = path_to_string(ghat, valid_paths[i])

					if path_lids[i] not in d_literal_paths:
						pe_idx += 1
						path_ids[i] = "PATHID_{:08d}".format(pe_idx)
						d_literal_paths[path_lids[i]] = path_ids[i]
						path_info = info_from_tids(dannot, dbundle['t'][i])
						num_paths += 1

						d_paths[path_ids[i]] = dict(
							gene_name=path_info[0], 
							gene_id=path_info[1],
							strand=path_info[2],
							tid=dbundle['t'][i], 
							ref=cid,
							exons=path_exon_list(ghat, valid_paths[i]), 
							in_event=set(), 
							in_loc=set())

					else:
						path_ids[i] = d_literal_paths[path_lids[i]]

					loc_path_set.append(path_ids[i])

				# bind paths together into events
				for asevent in bundle_results:

					aid = asevent[asevent[2]]
					bid = asevent[asevent[2]*-1+1]

					# first check for the literal event to be sure it hasn't been 
					# found already

					ae_lid = "{}-{}".format(path_lids[aid], path_lids[bid])
					if ae_lid not in d_literal_events:

						ae_idx += 1
						ae_id = "AEID_{:08d}".format(ae_idx)
						d_literal_events[ae_lid] = ae_id

						loc_ae_set.append(ae_id)
						ae_pair = [path_ids[aid], path_ids[bid]]

						d_events[ae_id] = dict(pid=ae_pair, type=asevent[3])
						num_events += 1

						# update the paths to note which ae event they are a part of
						d_paths[path_ids[aid]]['in_event'].update([ae_id])
						d_paths[path_ids[bid]]['in_event'].update([ae_id])

				# update loc. there may be multiple alt-events within a single locus. such as 
				# several alternative starts that all link into the same next exon. since 
				# everything is A/B then if we have 3 starts then we get 3 events all in 
				# one locus
				if len(loc_ae_set) > 0:
					aeloc_idx += 1
					aeloc_id = "AELOC_{:08d}".format(aeloc_idx)
					d_locs[aeloc_id] = dict(events=loc_ae_set, paths=loc_path_set)

					# update the paths
					for pid in loc_path_set:
						d_paths[pid]['in_loc'].update([aeloc_id])


	# final status for the bundle loop
	sys.stderr.write("processed {} bundles; {} paths in {} events in {} loci\n".format(tidx, num_paths, num_events, aeloc_idx))
	
	#
	# finished with transcript bundles
	#


	# this section deals with combining ae loci when paths appear in more than one locus.
	# this is accomplished by creating a graph
	sys.stderr.write("merging path clusters\n")
	
	g_loc = Graph(directed=False)
	# set a vertex up with a name attribute
#	g_loc.add_vertices(1)
#	g_loc.vs[0]['name'] = "foo"

	#
	# build edge list
	nlid = 0
	lida = lidb = 0
	ledges = [] # pairs of graph node indices connecting locid nodes
	dverts = {} # dict to keep track of the locid indices
	lverts = [] # list of the locid indices to apply as names for the graph nodes
	
	for pid in d_paths:
		# check if path appears in more than one locus
		llids = list(d_paths[pid]['in_loc'])
		nlids = len(llids)

		if nlids > 1:
			# need to combine these loci. check each one to see if it has been 
			# encountered already
			for locid in llids:
				if locid not in dverts:
					dverts[locid] = nlid
					lverts.append(locid)
					nlid += 1
			
			# create edges by just stringing these together. no need to make all 
			# combinations
			for i in range(1, nlids):
				ledges.append((dverts[llids[i-1]], dverts[llids[i]]))
	
	# add vertices and edges
	g_loc.add_vertices(nlid)
	g_loc.add_edges(ledges)
	# add name attribute to the vertices
	for i in range(nlid):
		g_loc.vs[i]['name'] = lverts[i]
			
	if False:
		# loop through path ids
		for pid in d_paths:
			# check if the path id appears in more than one alt-event locus. if so then
			# those loci must be merged.
			if len(d_paths[pid]['in_loc']) > 1:
				# get number of loci this appears in
				nlids = len(d_paths[pid]['in_loc'])
				# start list of locus ids. these will be used as nodes for the graph
				# we use to cluster all of these things
				lids = []
				# loop through locus ids
				for locid in list(d_paths[pid]['in_loc']):
					# see if the locus id already exists as a node in the graph
					rres = g_loc.vs.select(name_eq=locid)
					if len(rres)==0:
						g_loc.add_vertices(1)
						vid = g_loc.vcount()-1
						g_loc.vs[vid]['name'] = locid
					else:
						vid = rres[0].index
	
					lids.append(vid)
	
				# make pairwise connections
				for i in range(nlids-1):
					for j in range(1, nlids):
						g_loc.add_edges([(lids[i], lids[j])])

	# cluster the graph
	if g_loc.vcount() > 1:
		g_loc_clust = g_loc.clusters(mode=WEAK)
		g_loc_groups = g_loc_clust.membership
		g_loc_vnames = [v['name'] for v in g_loc.vs]
		dg_loc_groups = split(g_loc_vnames, g_loc_groups)

		# loop through the groups and merge the locs
		for loc_gid in dg_loc_groups.keys():
			lids = sorted(dg_loc_groups[loc_gid])
			if len(lids) > 1:
				# keep the first one, ditch the rest
				path_set = set(d_locs[lids[0]]['paths'])
				event_set = set(d_locs[lids[0]]['events'])

				for i in range(1, len(lids)):
					path_set.update(d_locs[lids[i]]['paths'])
					event_set.update(d_locs[lids[i]]['events'])
					d_locs[lids[i]] = None

				d_locs[lids[0]]['paths'] = list(path_set)
				d_locs[lids[0]]['events'] = list(event_set)
	

	if False:

#		for k in sorted(d_paths.keys()):
#			print k, d_paths[k]

#		for k in sorted(d_events.keys()):
#			print k, d_events[k]

		for k in sorted(d_locs.keys()):
			print k, d_locs[k]


	# I might export this in two ways one of them being the MISO format. 

	

	# each event in d_events is a "gene" and each path is an "mRNA" followed 
	# by each exon in the path


	# write an event table
	outfile = "{}.tsv".format(args.stub)
	sys.stderr.write("writing {}\n".format(outfile))
	fout = open(outfile, "w")

	fout.write("aeid\tae_type\tgene_name\tgene_id\ttranscript_id.A\ttranscript_id.B\tpath_id.A\tpath_id.B\tposition\n")

	for aeid in sorted(d_events.keys()):

		pidA = d_events[aeid]['pid'][0]
		pidB = d_events[aeid]['pid'][1]

		chrom = d_paths[pidA]['ref']
		strand = d_paths[pidA]['strand']
		estart = min([d_paths[pidA]['exons'][0][0], d_paths[pidB]['exons'][0][0]])
		eend = max([d_paths[pidA]['exons'][-1][1], d_paths[pidB]['exons'][-1][1]])

		gene_name = set(d_paths[pidA]['gene_name']).union(d_paths[pidB]['gene_name'])
		gene_id = set(d_paths[pidA]['gene_id']).union(d_paths[pidB]['gene_id'])

		lout = [aeid, 
			d_events[aeid]['type'], 
			",".join(list(gene_name)), 
			",".join(list(gene_id)), 
			",".join(list(d_paths[pidA]['tid'])),
			",".join(list(d_paths[pidB]['tid'])),
			d_events[aeid]['pid'][0], 
			d_events[aeid]['pid'][1], 
			"{}:{}-{}".format(chrom, estart, eend)]

		fout.write("\t".join(map(str, lout)) + "\n")

	fout.close()

	# --
	#
	# output the intron index
	#
	# --

	outfile = "{}_introns.tsv".format(args.stub)
	sys.stderr.write("writing {}\n".format(outfile))
	fout = open(outfile, "w")

	for pid in sorted(d_paths.keys()):
		nexons = len(d_paths[pid]['exons'])

		if nexons > 1:
			lints = []
			chrom = d_paths[pid]['ref']
			for i in range(1, nexons):
				lints.append("{}:{}-{}".format(chrom, 
					d_paths[pid]['exons'][i-1][1]+1, d_paths[pid]['exons'][i][0]-1))

			lout = [pid, ",".join(lints)]
			fout.write("\t".join(lout) + "\n")


	fout.close()



	# --
	#  
	# ouptut the GFF3
	#
	# --

	if args.miso:
		#
		# export MISO format
		#

		outfile = "{}.gff3".format(args.stub)
		sys.stderr.write("writing {}\n".format(outfile))
		gff_out = open(outfile, "w")		

		for k in sorted(d_events.keys()):
			ev = d_events[k]
			p1 = d_paths[ev['pid'][0]]
			p2 = d_paths[ev['pid'][1]]

			# need overall range
			gstart = p1['exons'][0][0]
			gend = p1['exons'][-1][1]

			if p2['exons'][0][0] < gstart:
				gstart = p2['exons'][0][0]

			if p2['exons'][-1][1] > gend:
				gend = p2['exons'][-1][1]

			gref = p1['ref']
			gstrand = p1['strand']
			ggid = list(set(p1['gene_id']).union(p2['gene_id']))
			ggname = list(set(p1['gene_name']).union(p2['gene_name']))

			id_gene = k
			id_mRna = ["{}.A".format(k), "{}.B".format(k)]


			# the gene row
			lout = [gref, 
				ev['type'], "gene", gstart, gend, ".", gstrand, ".", 
				"ID={};Name={};gid={}".format(id_gene, ",".join(ggname), ",".join(ggid))]

			#print "\t".join(map(str, lout))
			gff_out.write("\t".join(map(str, lout)) + "\n")

			# the first mRNA row

			pid = ev['pid'][0]
			lout = [gref, ev['type'], "mRNA", 
				p1['exons'][0][0], p1['exons'][-1][1], ".", gstrand, ".", 
				"ID={};Parent={};Name={};gid={}".format(id_mRna[0], id_gene, pid, ",".join(p1['gene_id']))]

			#print "\t".join(map(str, lout))
			gff_out.write("\t".join(map(str, lout)) + "\n")

			for i in range(len(p1['exons'])):
				lout = [gref, ev['type'], "exon", 
					p1['exons'][i][0], p1['exons'][i][1], ".", gstrand, ".", 
					"ID={}.{};Parent={};Name={}.exon{:02d};gid={}".format(id_mRna[0], i, 
						id_mRna[0], pid, i, ",".join(p1['gene_id']))]

				#print "\t".join(map(str, lout))
				gff_out.write("\t".join(map(str, lout)) + "\n")

			# the second mRNA row
			
			pid = ev['pid'][1]
			lout = [gref, ev['type'], "mRNA", 
				p2['exons'][0][0], p2['exons'][-1][1], ".", gstrand, ".", 
				"ID={};Parent={};Name={};gid={}".format(id_mRna[1], id_gene, pid, ",".join(p2['gene_id']))]

			#print "\t".join(map(str, lout))
			gff_out.write("\t".join(map(str, lout)) + "\n")

			for i in range(len(p2['exons'])):
				lout = [gref, ev['type'], "exon", 
					p2['exons'][i][0], p2['exons'][i][1], ".", gstrand, ".", 
					"ID={}.{};Parent={};Name={}.exon{:02d};gid={}".format(id_mRna[1], i, 
						id_mRna[1], pid, i, ",".join(p2['gene_id']))]

				#print "\t".join(map(str, lout))
				gff_out.write("\t".join(map(str, lout)) + "\n")

		gff_out.close()

	else:
		#
		# export GFF where each "gene" is an entire locus of paths. we'll have to 
		# export a second file that contains a map of paths to events.
		#

		outfile = "{}.gff3".format(args.stub)
		sys.stderr.write("writing {}\n".format(outfile))
		gff_out = open(outfile, "w")

		for aloc in d_locs.keys():

			if d_locs[aloc] is None:
				continue

			#
			# get start positions of all first exons from paths. also gather other
			# annotation from the paths
			#
			estarts = []
			eends = []
			gene_id = set()
			gene_name = set()
			strand = ""
			chrom = ""

			for pid in d_locs[aloc]['paths']:
				estarts.append(d_paths[pid]['exons'][0][0])
				eends.append(d_paths[pid]['exons'][-1][1])
				chrom = d_paths[pid]['ref']
				strand = d_paths[pid]['strand']
				gene_id.update(d_paths[pid]['gene_id'])
				gene_name.update(d_paths[pid]['gene_name'])

			gene_id = sorted(list(gene_id))
			gene_name = sorted(list(gene_name))
			pid_order = order(estarts)

			lout = [
				chrom, 
				".",
				"gene", 
				min(estarts), 
				max(eends), 
				".", strand, ".", 
				"ID={};Name={};Events={};".format(aloc, ",".join(gene_name), ",".join(d_locs[aloc]['events']))]

			#print "\t".join(map(str, lout))
			gff_out.write("\t".join(map(str, lout)) + "\n")

			# each path is an mRNA entry followed by each exon
			for i in pid_order:
				pid = d_locs[aloc]['paths'][i]
				
				lout = [chrom, ".", "mRNA", 
					d_paths[pid]['exons'][0][0], 
					d_paths[pid]['exons'][-1][1], 
					".", strand, ".", 
					"ID={};Parent={}".format(pid, aloc)]

				#print "\t".join(map(str, lout))
				gff_out.write("\t".join(map(str, lout)) + "\n")

				j = 0 # exon counter
				for exon in d_paths[pid]['exons']:
					j += 1

					lout = [chrom, ".", "exon", 
						exon[0], exon[1], 
						".", strand, ".", 
						"ID={}.exon{:02d};Parent={}".format(pid, j, pid)]

					#print "\t".join(map(str, lout))
					gff_out.write("\t".join(map(str, lout)) + "\n")

		gff_out.close()



#==============================================================================
# defs
#==============================================================================

def get_transcript_paths(annot, tid):

	ne = len(annot[tid])
	ledges = []

	i = 0
	last3 = ""
	for i in range(ne):
		# insert intron before next exon if we are past the first 
		# exon already
		if i > 0:
			ledges.append((last3, id5, "{}:{}".format(last3, id5), "i"))

		# start with the exon
		id3 = annot[tid][i].feature_3p_id() + ":3p"
		id5 = annot[tid][i].feature_5p_id() + ":5p"

		ledges.append((id5, id3, "{}:{}".format(id5, id3), "e"))


		last3 = id3

	return ledges


#
# returns sorted indices of values in v
def order(v):
	return np.argsort(v)

#
# get info for a set of transcript ids from the parsed GTF annotation
def info_from_tids(annot, tids):
	gnames = set()
	strand = ""
	gids = set()

	for tid in tids:
		gnames.update([annot[tid][0].attrs['gene_name']])
		gids.update([annot[tid][0].gid])
		strand = annot[tid][0].strand

	return (list(gnames), list(gids), strand)


# make a list of gene names from a list of transcript ids
def genes_from_tids(annot, tids):
	gnames = set()
	for tid in tids:
		gnames.update([annot[tid][0].attrs['gene_name']])

	return list(gnames)

def path_to_string(graph, path):
	vnames = [graph.vs[v]['name'] for v in path]
	return "@".join(vnames)

# make a list of exon start/end pairs from a given path within a given graph
def path_exon_list(graph, path):

	eset = []

	# first get any internal exons
	n = len(path)
	for i in range(1, n):
		eid = graph.get_eid(path[i-1], path[i])
		# get type
		if graph.es[eid]['type'] == 'e':
			tmp = map(int, graph.es[eid]['name'].split(":"))
			eset.append(tmp)

	if graph.vs[path[0]]['type'] == 'd':

		# first node of path is a donor so we need the upstream exon
		vups = graph.neighbors(path[0], mode="IN")
		vkeep = -1
		elen = 0
		for vid in vups:
			eid = graph.get_eid(vid, path[0])
			tmp = map(int, graph.es[eid]['name'].split(":"))
			tmplen = tmp[1]-tmp[0]+1
			if tmplen > elen:
				elen = tmplen
				vkeep = vid

		eid = graph.get_eid(vkeep, path[0])
		tmp = map(int, graph.es[eid]['name'].split(":"))

		eset = [tmp] + eset

	if graph.vs[path[-1]]['type'] == 'a':
		
		# last node is an acceptor so we need the next exon
		vdns = graph.neighbors(path[-1], mode="OUT")
		vkeep = -1
		elen = 0
		for vid in vdns:
			eid = graph.get_eid(path[-1], vid)
			tmp = map(int, graph.es[eid]['name'].split(":"))
			tmplen = tmp[1]-tmp[0]+1
			if tmplen > elen:
				elen = tmplen
				vkeep = vid

		eid = graph.get_eid(path[-1], vkeep)
		tmp = map(int, graph.es[eid]['name'].split(":"))

		eset.append(tmp)

	return eset

# 
# seta and setb are sets of vertices. this function will look at each pair
# and find all paths between them and validate the paths as it goes. for everything
# except alt starts and alt ends set min_tid to 2. for AS/AE then you need
# one common tid when finding valid paths from a start to an alt-in since
# that path will be uniquly from a single isoform, typically
def dig_for_paths(graph, seta, setb, min_tid=2):

	rres = []
	#
	# find paths, if any, between nodes in seta to nodes in setb. try all 
	# possible pairs. a vaoid path means all of the nodes in the path 
	# have at least one common transcript id
	for didx in seta:
		for aidx in setb:
			# do these nodes share transcript ids? (i.e. is it actually possible to have
			# valid paths from one to the other?)
			path_tid = graph.vs[didx]['tid'].intersection(graph.vs[aidx]['tid'])
			if len(path_tid) >= min_tid:

				# this means more than one transcript id is shared which suggests there 
				# may be more than one intron/exon path between the nodes
				lpaths = find_all_paths(graph, didx, aidx, path_tid)

				# make a set of valid paths from this list
				valid_paths = []
				valid_path_tid = []
				
				if min_tid==2:
				
					for path in lpaths:
						if len(path) > 0:
							final_path_tid, path_valid = validate_path(graph, path)
							if path_valid:
								valid_paths.append(path)
								valid_path_tid.append(final_path_tid)
	
					# do we have more than one valid path?
					if len(valid_paths) > 1:
						rres.append(dict(p=valid_paths, t=valid_path_tid))
						#gid_path_sets.append(dict(p=valid_paths, t=valid_path_tid))
				
				elif min_tid==1:
					# use this mode when looking for AS/AE events because each node-to-node
					# pair would probably only produce a single valid path. the calling code
					# will need to bundle the valid paths up but that's not the responsibility
					# of this function

					for path in lpaths:
						# only accept 3-node paths
						if len(path) == 3:
							final_path_tid, path_valid = validate_path(graph, path)
							if path_valid:
								valid_paths.append(path)
								valid_path_tid.append(final_path_tid)
	
					# do we have more than one valid path?
					if len(valid_paths) > 0:
						rres.append(dict(p=valid_paths, t=valid_path_tid))
						#gid_path_sets.append(dict(p=valid_paths, t=valid_path_tid))
					

	return rres

# 
# from the output of dig_for_paths when using it in AS/AE mode we can take the 
# list of valid paths and expand on it a bit. we have to take all unique end
# nodes of the path and then using those as a start point find all valid 
# 3-node paths that lead upstream (i.e. use the 3rd node as a alt-in and find
# incoming paths
def pair_AS(g, lpaths):

	unodes = set()
	upaths = set()
	npaths = []
	stmp = ""
	intronBundles = {}

	# create strings from the node sets that represent the ones that 
	# are from actual starts because when we pair them up we have to 
	# make sure that at least one of the pair is a start. 
	
	for i in range(len(lpaths)):
		# path to string
		stmp = ";".join(map(str, lpaths[i]['p'][0]))
		upaths.update([stmp])
		unodes.update([lpaths[i]['p'][0][2]])
	
	for nid in list(unodes):
		npaths = npaths + find_all_nstep_paths(g, nid, 3, mode="IN")
	
	# need to bundle these by the intron. since these are AS the first two nodes 
	# of the paths make up the a/d boundaries of an exon while the third is an a
	# boundary of another exon. so we need the d->a boundaries that are unique and
	# we want the one that has the longest exon
	for i in range(len(npaths)):
		p = npaths[i]['p'][0]
		stmp = ";".join(map(str, p[1:3]))
		utmp = ";".join(map(str, p))
		di = genome_dist_from_names(g.vs[p[0]]['name'], g.vs[p[1]]['name'])
		
		if stmp not in intronBundles:
			intronBundles[stmp] = {"len":0, "path":None, "flag":False}
			
			intronBundles[stmp]['len'] = di
			intronBundles[stmp]['path'] = npaths[i]
			
			if utmp in upaths:
				# this path is a start. we have to make sure we keep those
				# and do not replace them with another path that's not a start:
				# [=====]------[=====]-------[====] not a start
				#                 [==]-------[====] a start. don't lose these in this loop
				intronBundles[stmp]['flag'] = True
		else:
			# this one is setup already. we can update it if this is also a start or also 
			# not a start and it the exon is longer than the current on
			if di > intronBundles[stmp]['len']:
				if utmp in upaths:
					# this is start, update
					intronBundles[stmp]['len'] = di
					intronBundles[stmp]['path'] = npaths[i]
					intronBundles[stmp]['flag'] = True
				else:
					# this is not a start so only update if the current entry is also 
					# not a start
					if not intronBundles[stmp]['flag']:
						# update
						intronBundles[stmp]['len'] = di
						intronBundles[stmp]['path'] = npaths[i]
	
	#
	# final set of paths. these are the ones that we'll pair up to make AS events
	fpaths = []
	for pairId in intronBundles.keys():
		fpaths.append(dict(intronBundles[pairId]['path']))
	
	#
	# pair these up!
	n = len(fpaths)
	events = []
	for i in range(n-1):
		for j in range(i, n):
			p1 = fpaths[i]['p'][0]
			p2 = fpaths[j]['p'][0]
			
			if p1[2]==p2[2]:
				if p1[0] != p2[0] and p1[1] != p2[1]:
					# this is good!
					# make sure one of them is a start
					uid1 = ";".join(map(str, p1))
					uid2 = ";".join(map(str, p2))
					if (uid1 in upaths) or (uid2 in upaths):
						# good to go! set their order such that the one with the 
						# longer intron is first
						pos1 = genome_pos_from_name(g.vs[p1[1]]['name'])
						pos2 = genome_pos_from_name(g.vs[p2[1]]['name'])
						if pos1 < pos2:
							# first one is the first one
							events.append({'p':[fpaths[i]['p'][0], fpaths[j]['p'][0]], 't':[fpaths[i]['t'][0], fpaths[j]['t'][0]]})
						else:
							events.append({'p':[fpaths[j]['p'][0], fpaths[i]['p'][0]], 't':[fpaths[j]['t'][0], fpaths[i]['t'][0]]})

	return events


# 
# from the output of dig_for_paths when using it in AS/AE mode we can take the 
# list of valid paths and expand on it a bit. we have to take all unique end
# nodes of the path and then using those as a start point find all valid 
# 3-node paths that lead upstream (i.e. use the 3rd node as a alt-in and find
# incoming paths
def pair_AE(g, lpaths):

	unodes = set()
	upaths = set()
	npaths = []
	stmp = ""
	intronBundles = {}

	# create strings from the node sets that represent the ones that 
	# are from actual ends because when we pair them up we have to 
	# make sure that at least one of the pair is an end. 
	
	for i in range(len(lpaths)):
		# path to string
		stmp = ";".join(map(str, lpaths[i]['p'][0]))
		upaths.update([stmp])
		unodes.update([lpaths[i]['p'][0][0]])
	
	for nid in list(unodes):
		npaths = npaths + find_all_nstep_paths(g, nid, 3, mode="OUT")
	
	# need to bundle these by the intron. since these are AE the second two nodes 
	# of the paths make up the a/d boundaries of an exon while the first is a 'd'
	# boundary of another exon. so we need the d->a boundaries that are unique and
	# we want the one that has the longest exon
	for i in range(len(npaths)):
		p = npaths[i]['p'][0]
		stmp = ";".join(map(str, p[0:2]))
		utmp = ";".join(map(str, p))
		# get length of the exon
		di = genome_dist_from_names(g.vs[p[1]]['name'], g.vs[p[2]]['name'])
		
		if stmp not in intronBundles:
			intronBundles[stmp] = {"len":0, "path":None, "flag":False}
			
			intronBundles[stmp]['len'] = di
			intronBundles[stmp]['path'] = npaths[i]
			
			if utmp in upaths:
				# this path is a start. we have to make sure we keep those
				# and do not replace them with another path that's not a start:
				# [=====]------[=====]-------[====] not an end
				# [=====]------[==]                 an end. don't lose these in this loop
				intronBundles[stmp]['flag'] = True
		else:
			# this one is setup already. we can update it if this is also an end or also 
			# not an end and if the exon is longer than the current one
			if di > intronBundles[stmp]['len']:
				if utmp in upaths:
					# this is start, update
					intronBundles[stmp]['len'] = di
					intronBundles[stmp]['path'] = npaths[i]
					intronBundles[stmp]['flag'] = True
				else:
					# this is not a start so only update if the current entry is also 
					# not a start
					if not intronBundles[stmp]['flag']:
						# update
						intronBundles[stmp]['len'] = di
						intronBundles[stmp]['path'] = npaths[i]
	
	#
	# final set of paths. these are the ones that we'll pair up to make AE events
	fpaths = []
	for pairId in intronBundles.keys():
		fpaths.append(dict(intronBundles[pairId]['path']))
	
	#
	# pair these up!
	n = len(fpaths)
	events = []
	for i in range(n-1):
		for j in range(i, n):
			p1 = fpaths[i]['p'][0]
			p2 = fpaths[j]['p'][0]
			
			if p1[0]==p2[0]:
				# matching anchor
				if p1[1] != p2[1] and p1[2] != p2[2]:
					# this is good!
					# make sure one of them is an end
					uid1 = ";".join(map(str, p1))
					uid2 = ";".join(map(str, p2))
					if (uid1 in upaths) or (uid2 in upaths):
						# good to go! set their order such that the one with the 
						# longer intron is first
						pos1 = genome_pos_from_name(g.vs[p1[1]]['name'])
						pos2 = genome_pos_from_name(g.vs[p2[1]]['name'])
						if pos1 > pos2:
							# first one is the first one
							events.append({'p':[fpaths[i]['p'][0], fpaths[j]['p'][0]], 't':[fpaths[i]['t'][0], fpaths[j]['t'][0]]})
						else:
							events.append({'p':[fpaths[j]['p'][0], fpaths[i]['p'][0]], 't':[fpaths[j]['t'][0], fpaths[i]['t'][0]]})

	return events


#
# get position from node name
def genome_pos_from_name(n1):
	tmp = n1.split(":")
	p1 = float(tmp[0])
	return p1	

# 
# names of nodes are formatted as <pos>:<type>
def genome_dist_from_names(n1, n2):
	
	tmp = n1.split(":")
	p1 = float(tmp[0])
	tmp = n2.split(":")
	p2 = float(tmp[0])
	return p2-p1+1

#
# compares a pair of paths and figures out if the pair represent one of the 
# typical alternative event types
def compare_paths(graph, a, b):

	# variables
	rres = None
	ref = -1

	#
	# get lengths
	lenA = len(a)
	lenB = len(b)
	#
	# get start and end types
	start_type = graph.vs[a[0]]['type']
	end_type = graph.vs[a[-1]]['type']
	#
	# check if starts and ends are the same (i.e. we have a box). 'box' types
	# include cassettes, mutually ex, alt 3', alt 5' and retained intron.
	same_start = a[0]==b[0]
	same_end = a[-1]==b[-1]

	# first check lengths.  if they are both longer than 4 then this isn't a 
	# simple alt path
	if lenA > 4 and lenB > 4:
		return (rres, ref)

	# lengths are OK, now check if we have a box type
	if same_start and same_end:

		if start_type=='d' and end_type=='a':
			#
			# box type from alt donor to alt acceptor. 
			# check lengths to determine what is going on. for SE and TE we 
			# call the path with the included exon(s) the reference so the resulting 
			# quantification is percent spliced in (PSI) as the convention goes
			if lenA == 2 and lenB > 2:
				ref = 1
				# skip type
				if lenB == 4:
					rres = "SE"
				elif lenB > 4:
					rres = "TE"

			elif lenA > 2 and lenB == 2:
				ref = 0
				if lenA == 4:
					rres = "SE"
				elif lenA > 4:
					rres = "TE"

			elif lenA == 4 and lenB == 4:
				# possibly mutually exclusive. confirm that the exon is not shared in any 
				# transcript
				eidA = graph.get_eid(a[1], a[2])
				eidB = graph.get_eid(b[1], b[2])

				tid_int = graph.es[eidA]['tid'].intersection(graph.es[eidB]['tid'])
				if len(tid_int)==0:
					# good! now make sure the exons do not overlap at all
					
					# get exon boundaries
					boundA = map(float, graph.es[eidA]['name'].split(":"))
					boundB = map(float, graph.es[eidB]['name'].split(":"))

					if not boundary_overlap(boundA, boundB):
						# this is a mutually exlusive event
						rres = "ME"
						# always use the path with the ME exon that's closest to the 
						# alt donor site
						ref = 0
						if boundB[0] < boundA[0]:
							ref = 1

		elif start_type=='a' and end_type=='a':
			# acceptor to acceptor is alt 5' so both paths have to be 3 nodes
			if lenA==3 and lenB==3:
				# good!
				rres = "A5"
				# ref should be the one with the longer intron (or shorter exon)
				eidA = graph.get_eid(a[0], a[1])
				eidB = graph.get_eid(b[0], b[1])
				tmp = map(int, graph.es[eidA]['name'].split(":"))
				elenA = tmp[1]-tmp[0]+1
				tmp = map(int, graph.es[eidB]['name'].split(":"))
				elenB = tmp[1]-tmp[0]+1

				if elenA > elenB:
					ref = 1
				else:
					ref = 0 

		elif start_type=='d' and end_type=='d':
			# donor to donor is alt 3' so both paths have to be 3 nodes
			if lenA==3 and lenB==3:
				# good!
				rres = "A3"
				# ref should be the one with the longer intron (or shorter exon)
				eidA = graph.get_eid(a[1], a[2])
				eidB = graph.get_eid(b[1], b[2])
				tmp = map(int, graph.es[eidA]['name'].split(":"))
				elenA = tmp[1]-tmp[0]+1
				tmp = map(int, graph.es[eidB]['name'].split(":"))
				elenB = tmp[1]-tmp[0]+1

				if elenA > elenB:
					ref = 1
				else:
					ref = 0 

		elif start_type=='a' and end_type=='d':
			# acceptor to donor - retained intron type. one path should be 2 nodes 
			# and the other should be 4. this is like an SE event but it happens 
			# within an exon instead of within an intron. reference is always the 
			# path without the intron so the psi value is percent retained.
			if lenA==2 and lenB==4:
				rres = "RI"
				ref = 0
			elif lenA==4 and lenB==2:
				rres = "RI"
				ref = 1

	elif same_start and not same_end:
		# alt end type with common start node but different end nodes. these paths
		# would only be 3 nodes long
		if lenA==3 and lenB==3:
			# one or both of these should be an actual end and also the center point
			# should not be the same (i.e. then they'd have the same intron)
			if a[1] != b[1]:
				# get outgoing list for a[2] and b[2]. one of the two needs to have 
				# no outgoing connections
				a_out = graph.neighbors(a[2], mode="OUT")
				b_out = graph.neighbors(b[2], mode="OUT")

				if len(a_out) == 0 or len(b_out)==0:
					# good to go
					rres = "AE"
					if len(a_out)==0:
						ref = 0
					else:
						ref = 1

	elif same_end and not same_start:
		# alt start type with common end node but different start nodes. these paths
		# would only be 3 nodes long
		if lenA==3 and lenB==3:
			# one or both of these should be an actual start and also the center point
			# should not be the same (i.e. then they'd have the same intron)
			if a[1] != b[1]:
				# get incoming list for a[0] and b[0]. one of the two needs to have 
				# no incoming connections
				a_in = graph.neighbors(a[0], mode="IN")
				b_in = graph.neighbors(b[0], mode="IN")

				if len(a_in) == 0 or len(b_in)==0:
					# good to go
					rres = "AS"
					if len(a_in)==0:
						ref = 0
					else:
						ref = 1

	return (rres, ref)

def boundary_overlap(a, b):
	rres = False
	if a[1] > b[0] and a[0] < b[1]:
		rres = True

	return rres

#
# this function validates paths. for a path to be valid the nodes must alternate
# between 'd' and 'a' types, starting from the first one, and the intersection
# of all of the tid's from all nodes and edges must leave at least 1 transcript
def validate_path(graph, path):
	rres = True

	# loop through the path
	if len(path) < 2:
		return [], False

	ltype = graph.vs[path[0]]['type']
	tidset = graph.vs[path[0]]['tid']

	for i in range(1, len(path)):
		ctype = graph.vs[path[i]]['type']
		if ctype == ltype:
			rres = False
			break
		ltype = ctype

		# get edge id
		eid = graph.get_eid(path[i-1], path[i])
		# check tid set for the edge
		tidset0 = tidset.intersection(graph.es[eid]['tid'])
		if len(tidset0) < 1:
			rres = False
			break

		tidset = tidset0

		# check tid intersection
		tidset0 = tidset.intersection(graph.vs[path[i]]['tid'])
		if len(tidset0) < 1:
			rres = False
			break

		tidset = tidset0

	return tidset, rres

def find_all_paths(graph, source, target, path_tid, mode="OUT"):
	# get adj node list
	adjout = graph.get_adjlist(mode=mode)
	vtid = [v['tid'] for v in graph.vs]
	return adjlist_find_paths(adjout, vtid, source, target, path_tid)

def adjlist_find_paths(a, vtid, s, t, ptid, path=[]):
	# append current "from" node to path if it has intersection with the ptid set
	pint = ptid.intersection(vtid[s])

	if len(pint) > 0:
		path = path + [s]
	else:
		# this path is done because we don't have any transcript id overlap with 
		# those that should be included
		return [[]]

	if s==t:
		# done, made it to t
		return [path]
	elif len(path) > (MAX_SKIPS*2 + 2):
		# path is getting too long - don't care about massive numbers of skipped
		# exons
		return [[]]

	if len(a[s])==0:
		# this one didn't make it
		return [[]]

	# continue the search!
	paths = []
	for child in a[s]:
		if child not in path:
			child_paths = adjlist_find_paths(a, vtid, child, t, ptid, path)
			for child_path in child_paths:
				paths.append(child_path)

	return paths


def find_all_nstep_paths(graph, source, nsteps, mode="OUT"):
	# get adj node list
	adjout = graph.get_adjlist(mode=mode)
	vtid = [v['tid'] for v in graph.vs]
	# initial tid set
	path_tid = graph.vs[source]['tid']
	#return adjlist_find_nstep_paths(adjout, vtid, source, nsteps, path_tid)
	rres0 = adjlist_find_nstep_paths_dev(adjout, source, nsteps, [])
	rres = []
	for p in rres0:
		
		if mode=="IN":
			# in "IN" mode the paths need to be reversed so they represent
			# the usual direction of the transcripts
			p.reverse()

		tset, valid = validate_path(graph, p)
		if valid:
			rres.append({"p":[p], "t":[tset]})
			
	return rres

def adjlist_find_nstep_paths(a, vtid, s, n, ptid, path=[]):
	# append current "from" node to path if it has intersection with the ptid set
	pint = ptid.intersection(vtid[s])

	if len(pint) > 0:
		path = path + [s]
	else:
		# this path is done because we don't have any transcript id overlap with 
		# those that should be included
		return [[]]

	if len(path) == (n+1):
		# done, made it to number of steps
		return [path]

	if len(a[s])==0:
		# this one didn't make it
		return [[]]

	# continue the search!
	paths = []
	for child in a[s]:
		if child not in path:
			child_paths = adjlist_find_nstep_paths(a, vtid, child, n, ptid, path)
			for child_path in child_paths:
				paths.append(child_path)

	return paths

def adjlist_find_nstep_paths_dev(alist, vid, n, p0):
	
	p0 = p0 + [vid]
	
	if len(p0)==n:
		# done
		return [p0]
	
	if len(alist[vid]) == 0:
		# done
		return [p0]
	
	p = []
	for cid in alist[vid]:
		if cid not in p0:
			ptmp = adjlist_find_nstep_paths_dev(alist, cid, n, p0)
			for ctmp in ptmp:
				p.append(ctmp)
		
	
	return p
	

def split(subject, factor):
	# loop through factor and split into a dict
	dout = {}
	for i in range(len(factor)):
		fid = "{}".format(factor[i])
		if fid not in dout:
			dout[fid] = []
		dout[fid].append(subject[i])

	return dout

def exon_id_set(grow):
	p5 = "{}:A".format(grow.start)
	p3 = "{}:D".format(grow.end)
	eid = "{}:{}".format(grow.start, grow.end)
	return (p5, p3, eid)

def intron_id_set(g1, g2):
	p5 = "{}:D".format(g1.end)
	p3 = "{}:A".format(g2.start)
	iid = "{}:{}".format(g1.end, g2.start)
	return (p5, p3, iid)

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
		return sz

	def feature_3p_id(self):
		sz = "{}:{}".format(self.rname, self.end)
		return sz


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Build alternative splicing index GFF (like MISO style)")
parser.add_argument('gtf', type=str, 
	help="GTF to base index on.")
parser.add_argument('stub', type=str, 
	help="Stub for output files")
parser.add_argument('--miso', action="store_const", const=True, default=False, 
	help="Output GFF in MISO format")
parser.add_argument('--introns', action="store_const", const=True, default=False, 
	help="Output a per-path intron index so that the events may be quantified from junctions")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

