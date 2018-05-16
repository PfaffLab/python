#!/usr/bin/env python
#==============================================================================
# seal-mappability-reduce.py
#
# Shawn Driscoll
# 6 March 2017
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script is meant to run on the output of seal-mappability-graph.py. 
# We're going to try to collect some useful statistics about the mapability 
# info and make it functionally useful in terms of a gene annotation (GTF) 
# which must be provided. Match scores between isoforms of a gene will be 
# considered differently than match scores between features that do not 
# share gene ids. This way we can learn something at the gene id level
# about how 'mappable' a gene id is relative to other gene ids. 
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib
from math import log,exp

# from subprocess import Popen
# from random import gauss, random, sample
import scipy.stats as scs
import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

#==============================================================================
# main
#==============================================================================

def main(args):

	hpares = {}
	lnomatch = []
	kswitch = {} # this is going to be a translation of the pair names
	kgswitch = {}
	khash = {}
	inter_gid = {}
	intra_gid = {}
	intra_gid_sim = {}
	intra_scores = []
	inter_scores = []
	tid2kid = {}
	tidmpp = {}
	gidmpp = {}
	tid2gid = {}
	gid2tid = {}
	htidraw = {}

	tidGraph = {}
	gidGraph = {}

	gid2kid = {}

	klen = None
	dgtf = None
	dfai = {}
	tid2k = {}

	# parse the entire mappability thing into a hash
	hmapp = {}

	finname = args.file
	tmp = finname.split("/")
	out_stub = tmp[-1]

	#
	# open and load the GTF
	try:
		message("Loading annotation from {}".format(args.gtf))
		dgtf = parse_gtf(args.gtf)
		message("Loaded {} transcripts".format(len(dgtf.keys())))
	except:
		return 1

	#
	# open and load the FAI
	message("Loading transcript lengths from {}".format(args.fai))
	fin = open(args.fai, "r")
	
	for szl in fin:
		aln = szl.strip().split("\t")
		dfai[aln[0]] = float(aln[1])

	fin.close()

	#
	# build tables from the GTF for converting tids to gids, etc
	for tid in dgtf.keys():
		gid = dgtf[tid][0].gid

		if tid not in tid2gid:
			tid2gid[tid] = gid

		if gid not in gid2tid:
			gid2tid[gid] = set()

		gid2tid[gid].update([tid])


	message("Loading mapability from {}".format(args.file))
	fin = open(args.file, "r")
	# skip header
	tmp = fin.readline()

	# read this thing in and split the connections into some 
	# buckets
	for szl in fin:
		aln = szl.strip().split("\t")

		tid0 = aln[0]
		tid1 = aln[1]
		shared_ratio = float(aln[2])
		shared_k_num = float(aln[3])
		ref_len = float(aln[4])
		if klen is None:
			klen = float(aln[5])
		ref_num_k = float(aln[6])

		# keep the number of kmers generated from this transcript
		if tid0 not in tid2k:
			tid2k[tid0] = ref_num_k

		if tid1 == "none":
			lnomatch.append(szl)
			continue

		# keep the raw number of shared kmers
		if tid0 not in hmapp:
			hmapp[tid0] = {}

		hmapp[tid0][tid1] = shared_k_num


#		k = make_keys(aln)
#		if k[0] not in kswitch:
#			# make it so either of the possible pairs both point to the 
#			# same thing
#			kswitch[k[0]] = k[0]
#			kswitch[k[1]] = k[0]
#			
#		kid = kswitch[k[0]]
#		
#		if kid not in khash:
#			khash[kid] = []
#		
#		khash[kid].append(list(aln))
#		
#		if tid0 not in tid2kid:
#			tid2kid[tid0] = set()
#		
#		tid2kid[tid0].update([kid])
		
		if tid1 != "none":
			#
			# gene gene ids for the two transcripts
			#
			gid1 = tid2gid[tid0]
			gid2 = tid2gid[tid1]

			#
			# initalize entry in the gid graph and tid graph for this connection
			#

			#
			# establish key for the distinct trainscript pair
			# 

			k = make_keys([tid0, tid1])

			if (k[0] in kswitch) and (k[1] in kswitch):
				ktid = kswitch[k[0]]
			else:
				# not both in there
				if (k[0] in kswitch) or (k[1] in kswitch):
					# why would only one of them be present?
					sys.stderr.write("wtf, bro\n")
					if k[0] in kswitch:
						ktid = kswitch[k[0]]
						kswitch[k[1]] = ktid
					else:
						ktid = kswitch[k[1]]
						kswitch[k[0]] = ktid
				else:
					# neither are present
					ktid = k[0]
					kswitch[k[0]] = ktid
					kswitch[k[1]] = ktid

			if ktid not in tidGraph:
				# make entry for this pair. we will need the sim scores and the feature lengths
				tidGraph[ktid] = {'sim':aln[3], 'len':[dfai[tid0]]}
			else:
				tidGraph[ktid]['sim'].append(aln[2])
				tidGraph[ktid]['len'].append(dfai[tid0])
				print ktid, tidGraph[ktid]


			# 
			# initalize an entry for this transcript
			if tid0 not in tidmpp:
				tidmpp[tid0] = {'intra':[], 'inter':[], 'intra_tid':[], 'inter_tid':[]}

			#
			# is this a connection that's inter-gid (between) or intra-gid (within)?
			#

			if gid1==gid2:
				
				# intra means the connection is between transcripts of the same 
				# gene id. we'll monitor this in terms of gene ids

				if gid1 not in intra_gid:
					intra_gid[gid1] = set()
				
				kid = "{}|{}".format(tid0, tid1)

				intra_gid[gid1].update([kid])

				if gid1 not in intra_gid_sim:
					intra_gid_sim[gid1] = []

				# insert the sim score to this gene id
				intra_gid_sim[gid1].append(aln[2])

				# collect the intra connection scores
				intra_scores.append(aln[2])
				tidmpp[tid0]['intra'].append(aln[2])
				tidmpp[tid0]['intra_tid'].append(tid1)
				
			else:

				# inter means the connection is between transcripts of different 
				# gene ids. we'll monitor this in terms of the connection pair
				kid = "{}|{}".format(tid0, tid1)

				tidmpp[tid0]['inter'].append(aln[2])
				tidmpp[tid0]['inter_tid'].append(tid1)

				kg = make_keys([gid1, gid2])
				if kg[0] not in kgswitch:
					# make both forms point to the same id
					kgswitch[kg[0]] = kg[0]
					kgswitch[kg[1]] = kg[0]

				kgid = kgswitch[kg[0]]

				if kgid not in inter_gid:
					inter_gid[kgid] = set()

				inter_gid[kgid].update([kid])

	fin.close()

	#
	# summarize the intra-gene id transcript to transcript connections 
	# using the 'pair_sim' score
	message("Calculating typical intra-gid similarity")
	intra_scores = []
	pdone = set()
	withinGidSim = {}
	for gid in gid2tid.keys():
		
		withinGidSim[gid] = []

		if len(gid2tid[gid]) > 1:
			# multi transcript feature
			tset = list(gid2tid[gid])
			# loop through transcripts in this gene id and see if they have scores with their
			# intra gene id brothers
			for i in range(len(tset)-1):
				tid0 = tset[i]
				for j in range(i+1, len(tset)):
					tid1 = tset[j]
					k = make_keys([tid0, tid1])
					# skip if we already have this pair
					if (k[0] in pdone) or (k[1] in pdone):
						continue

					if tid0 in hmapp and tid1 in hmapp:
						if tid1 in hmapp[tid0] and tid0 in hmapp[tid1]:
							# add the pair keys into pdone
							pdone.update(k)
							# both connections are present
							s1 = hmapp[tid0][tid1]
							s2 = hmapp[tid1][tid0]
							l1 = dfai[tid0]
							l2 = dfai[tid1]
							v = pair_sim(s1, s2, l1, l2)
							intra_scores.append(v)

							# track all wihtin gid sims
							withinGidSim[gid].append(v)

	#
	# summarize inter-gid tid to tid similiarites using 'pair_sim'
	message("Calculating typical inter-gid similarity")
	inter_scores = []
	pdone = set()
	
	# hash to track all gene id to gene id connections by id
	interGidSim = {}

	# graph with keys for each distinct pair of gene ids
	ggraph = {}

	for tid in tidmpp.keys():
		if len(tidmpp[tid]['inter']) > 0:

			# this transcript has inter-gid connections
			
			gid0 = tid2gid[tid]
			v = []
			
			for i in range(len(tidmpp[tid]['inter'])):
				tid1 = tidmpp[tid]['inter_tid'][i]
				# make keys
				k = make_keys([tid, tid1])
				# has this pair been recorded yet?
				if (k[0] in pdone) or (k[1] in pdone):
					continue

				gid1 = tid2gid[tid1]
				if gid0 != gid1:
					# inter connection. check to see if the opposite connection is present...
					if tid1 in hmapp:
						if tid in hmapp[tid1]:
							# add the pair to pdone
							pdone.update(k)

							# yes, get the two scores
							s1 = hmapp[tid][tid1]
							s2 = hmapp[tid1][tid]
							l1 = dfai[tid]
							l2 = dfai[tid1]
							v = pair_sim(s1, s2, l1, l2)

							inter_scores.append(v)

							if gid0 not in interGidSim:
								interGidSim[gid0] = {}

							if gid1 not in interGidSim[gid0]:
								interGidSim[gid0][gid1] = []

							interGidSim[gid0][gid1].append(v)


	#
	# Global level stats on the inter- and intra-gid shared similarites. mappatilibites would be 
	# 1-x where 'x' is either of the means. the variances would be identical. 
	# mean and sd of the intra connection scores
	intra_mean = np.mean(intra_scores)
	intra_std = np.std(intra_scores)
	inter_mean = np.mean(inter_scores)
	inter_std = np.std(inter_scores)

	#
	# in this loop build a new hash that's keyed by gid and at each entry we have a list
	# of all distinct sim scores associated with it. so just the gid-2-gid averages 
	# that are generated in this loop. then we'll use the average of the average
	#
	

	pdone = set()
	message("Creating gid-gid graph")

	fout = open("{}.gid_graph".format(out_stub), "w")
	fout.write("\t".join(["from", "to", "shared", "inter_rank", "intra_rank", "inter_prob", "intra_prob", "inter_intra_llr"]))
	fout.write("\n")

	for gid0 in gid2tid.keys():
		
		if gid0 in interGidSim:
			for gid1 in interGidSim[gid0].keys():
				
				k = make_keys([gid0, gid1])
				if (k[0] in pdone) or (k[1] in pdone):
					continue

				pdone.update(k)

				# get values
				v = interGidSim[gid0][gid1]
				# check for the opposite connection
				if gid1 in interGidSim:
					if gid0 in interGidSim[gid1]:
						v += interGidSim[gid1][gid0]

				mu = np.mean(v)
				vv = max([np.std(v)**2, inter_std**2])

				inter_rank = (mu-inter_mean)/math.sqrt(inter_std**2 + vv)
				inter_prob = scs.norm.cdf(inter_rank, 0, 1)
				intra_rank = (mu-intra_mean)/math.sqrt(intra_std**2 + vv)
				intra_prob = scs.norm.cdf(intra_rank, 0, 1)
				interIntraLLR = log(inter_prob/intra_prob)

				l0 = [gid0, gid1]
				lvals = ["{:0.6f}".format(x) for x in [mu, inter_rank, intra_rank, inter_prob, intra_prob, interIntraLLR]]
				
				fout.write("\t".join(l0+lvals))
				fout.write("\n")

	fout.close()



	return 1

#	print len(intra_scores), intra_mean, intra_std, len(inter_scores), inter_mean, inter_std
#	return 1

	intra_mean_var = []
	inter_mean_var = []

	message("Calculating typical intra- and inter-gid mappability for each transcript")
	for tid in tidmpp:
		#
		# stats
		intra_stats = None
		inter_stats = None

		if len(tidmpp[tid]['intra']) > 0:
			# has intra. take maximum shared value
			v = [1-x for x in tidmpp[tid]['intra']]
			intra_stats = { 'max':max(v), 'min':min(v), 'mean':np.mean(v), 'stdev':np.std(v) }
			if intra_stats['stdev'] > 0:
				intra_mean_var.append(intra_stats['stdev']**2)

		if len(tidmpp[tid]['inter']) > 0:
			# has inter connections. take maximum shared value and typical shared value.
			v = [1-x for x in tidmpp[tid]['inter']]
			inter_stats = { 'max':max(v), 'min':min(v), 'mean':np.mean(v), 'stdev':np.std(v) }
			if inter_stats['stdev'] > 0:
				inter_mean_var.append(inter_stats['stdev']**2)

		tidmpp[tid]['intra_stats'] = intra_stats
		tidmpp[tid]['inter_stats'] = inter_stats

#		if (intra_stats is not None) and (inter_stats is not None):
#			foo = (inter_stats['mean']-intra_stats['mean'])/math.sqrt(inter_stats['stdev']**2 + intra_stats['stdev']**2)
#			print foo

	
	# denom sqrt for comparing intra- to inter-gid values
	tid2tid_intra_mapp_var = np.mean(intra_mean_var)
	tid2tid_inter_mapp_var = np.mean(inter_mean_var)

	for tid in tidmpp:

		tid_num = intra_mean
		if (tidmpp[tid]['intra_stats'] is not None):
			tid_num = tidmpp[tid]['intra_stats']['mean']

		if (tidmpp[tid]['intra_stats'] is not None) and (tidmpp[tid]['inter_stats'] is not None):
			# calculate inter to intra score
			foo = (tidmpp[tid]['inter_stats']['mean']-tidmpp[tid]['intra_stats']['mean'])/foo_denom
			print tid, foo


	#
	# assemble information about gene id to gene id connections
	#

	message("Building gene id to gene id connection data")

	# gid to gid hash
	hgid2gid = {}

	for gid in gid2tid.keys():
		
		tset = list(gid2tid[gid])
		conn = False
		if gid not in hgid2gid:
			hgid2gid[gid] = {}
		
		# loop through transcripts in this gid to find connections
		for tid in tset:
			# check if the tid has connections
			if tid in tidmpp:
				# does it have 'inter' connections?
				if len(tidmpp[tid]['inter']) > 0:
					# yes
					conn = True
					# loop through transcript connections and dump them out 
					# into gid buckets. these are 'directional' connections because
					# we want to track the gene id's mappability relative to 
					# other gene ids.
					for i in range(len(tidmpp[tid]['inter'])):
						ctid = tidmpp[tid]['inter_tid'][i]
						cgid = tid2gid[ctid]
						csim = tidmpp[tid]['inter'][i]

						if cgid not in hgid2gid[gid]:
							hgid2gid[gid][cgid] = {'weights':[], 'sims':[] }

						# append the score
						hgid2gid[gid][cgid]['sims'].append(csim)
						# append the length of the current 'tid' as the weight
						hgid2gid[gid][cgid]['weights'].append(dfai[tid])

	#
	# now we should have a hash table of all gid keys and some of them may have 
	# connections to other gids. those connections will have a set of 
	# values from each connection recorded.
	#

	message("Summarizing into an edge table and writing to file")

	fout = open("{}.gid_graph".format(out_stub), "w")
	lgid2gid = []
	for gid in hgid2gid.keys():
		kgid = hgid2gid[gid].keys()
		
		ltmp = [gid, "none", "NA", "NA"]

		if len(kgid) > 0:
			# we have connection(s)
			# loop through them
			for cgid in kgid:
				# create a summary
				ltmp[1] = cgid
				# get mean sim. the score are weighted by the length of the transcript of this gid
				# that the connection was based on.
				ltmp[2] = wmean(hgid2gid[gid][cgid]['sims'], hgid2gid[gid][cgid]['weights'])
				#ltmp[3] = log(scs.norm.cdf(ltmp[2], intra_mean, intra_std)/0.5)
				ltmp[3] = (ltmp[2]-intra_mean)/intra_std
				#print gid, cgid, ltmp[2], hgid2gid[gid][cgid]['sims'], hgid2gid[gid][cgid]['weights']

				#lgid2gid.append(list(ltmp))
				fout.write("\t".join(map(str, ltmp)))
				fout.write("\n")
		else:
			fout.write("\t".join(map(str, ltmp)))
			fout.write("\n")

	fout.close()

	#
	# repeat the code from above that generated tidmpp to make gidmpp
	#

	message("Building per gene id summary for those that have connections")

	gidNoMatch = []
	gidmpp = {}
	fin = open("{}.gid_graph".format(args.file), "r")

	for szl in fin:
		aln = szl.strip().split("\t")

		if aln[1] == "none":
			gidNoMatch.append("\t".join(map(str, aln)))
			continue

		aln[2] = float(aln[2])
		aln[3] = float(aln[3])

		# this is a between gid connection

		if aln[0] not in gidmpp:
			gidmpp[aln[0]] = {'inter':[], 'inter_gid':[]}
			
		gidmpp[aln[0]]['inter'].append(aln[2])
		gidmpp[aln[0]]['inter_gid'].append(aln[1])

	fin.close()

	#
	# now we can generate a gene_id level report
	#

	message("Writing gene_id summary table to {}.gid_summary".format(out_stub))
	fout = open("{}.gid_summary".format(out_stub), "w")

	for gid in sorted(gid2tid.keys()):

		minter = "1"
		minter_prob = "NA"
		inter_gid = "None"
		inter_gid_vals = "None"
		lout = []

		if gid in gidmpp:

			tmp = geom_mean(flip_v(gidmpp[gid]['inter']))
			minter_prob = "{:f}".format(log(scs.norm.cdf(1-tmp, intra_mean, intra_std)/0.5))
			minter = "{:0.4f}".format(tmp)
			inter_gid = ",".join(gidmpp[gid]['inter_gid'])

			inter_gid_vals = ""
			for i in range(len(gidmpp[gid]['inter'])):
				inter_gid_vals += "{:0.4f},".format(gidmpp[gid]['inter'][i])


		lout = [gid, "{:d}".format(len(gid2tid[gid])), 
			minter, 
			minter_prob,
			inter_gid, 
			re.sub(",$", "", inter_gid_vals)]

		fout.write("\t".join(lout) + "\n")

	fout.close()


	#
	# write output
	message("Writing transcript_id summary table to {}.tid_summary".format(out_stub))
	fout = open("{}.tid_summary".format(out_stub), "w")

	#
	# this section prints out a transcript mapability summary
	for tid in sorted(dgtf.keys()):
		minter = "NA"
		mintra = "NA"
		minter_prob = "NA"
		mintra_prob = "NA"
		inter_tid = "None"
		inter_tid_vals = "None"
		intra_tid = "None"
		intra_tid_vals = "None"
		lout = []
		v0 = []

		#
		# get the gene id for this transcript and then make a hash of the 
		# other transcripts wihin it's gene id, if any
		gid = tid2gid[tid]
		stid = gid2tid[gid]

		# this section deals with 'intra' only. 'inter' is handled next
		if len(stid) > 1:
			# multiple
			htid = {}
			# initalize hash of the tids with what will be their 
			# shared ratios if they are extracted below
			for k in list(stid):
				if k != tid:
					htid[k] = "NA"

			if tid in tidmpp:
				if len(tidmpp[tid]['intra']) > 0:
					# get summary
					tmp = geom_mean(flip_v(tidmpp[tid]['intra']))
					mintra_prob = "{:f}".format(log(scs.norm.cdf(1-tmp, intra_mean, intra_std)/0.5))
					mintra = "{:0.4f}".format(tmp)
					# we have intra-gid connections for this transcript
					for i in range(len(tidmpp[tid]['intra'])):
						# copy values
						htid[tidmpp[tid]['intra_tid'][i]] = "{:0.4f}".format(tidmpp[tid]['intra'][i])

					# build the strings to print
					ktid = htid.keys()
					intra_tid = ktid[0]
					intra_tid_vals = htid[ktid[0]]
					if len(ktid) > 1:
						for i in range(1, len(ktid)):
							intra_tid += ",{}".format(ktid[i])
							intra_tid_vals += ",{}".format(htid[ktid[i]])

				else:
					# compared to inter gid but not intra
					mintra = "NA"
					intra_tid_vals = ",".join(["NA" for i in range(len(stid)-1)])
					intra_tid = ",".join(htid.keys())

			else:
				# not compared to anything or at least not reported from 'seal-mapability-graph'
				mintra = "NA"
				intra_tid_vals = ",".join(["NA" for i in range(len(stid)-1)])
				intra_tid = ",".join(htid.keys())


		else:
			# not more than 1 transcript in this gene id
			intra_tid = "None"
			intra_tid_vals = "None"
			mintra = "NA"



		if tid in tidmpp:

			if len(tidmpp[tid]['inter']) > 0:
				tmp = geom_mean(flip_v(tidmpp[tid]['inter']))
				minter_prob = "{:f}".format(log(scs.norm.cdf(1-tmp, intra_mean, intra_std)/0.5))
				minter = "{:0.4f}".format(tmp)
				inter_tid = ",".join(tidmpp[tid]['inter_tid'])

				inter_tid_vals = ""
				for i in range(len(tidmpp[tid]['inter'])):
					inter_tid_vals += "{:0.4f},".format(tidmpp[tid]['inter'][i])


		lout = [tid, tid2gid[tid], "{:d}".format(len(stid)), 
			mintra, 
			minter, 
			mintra_prob,
			minter_prob,
			intra_tid,
			inter_tid, 
			re.sub(",$", "", intra_tid_vals), 
			re.sub(",$", "", inter_tid_vals)]

		fout.write("\t".join(lout) + "\n")

	fout.close()
	
	return 0

def make_keys(v):
	k1 = "{}|{}".format(v[0], v[1])
	k2 = "{}|{}".format(v[1], v[0])

	return(list([k1, k2]))

# collect the shared ratios from each pair in the set and generate a summary value
def joined_summary(l, mode=0):
	n = len(l)
	v = []
	rres = 0

	if n == 1:
		return(l[0][2])
	else:
		for i in range(n):
			v.append(l[i][2])

	if mode == 0:
		# return mean
		rres = mean(v)
	elif mode==1:
		# harmonic mean
		rres = harmonic_mean(v)

	return rres

def flip_v(v):
	n = len(v)
	vhat = []

	for i in range(n):
		vhat.append(1-v[i])

	return vhat

#
# weighted mean
def wmean(v, w):

	n = len(v)
	if n != len(w):
		return 0

	nsum = 0
	ndenom = 0
	for i in range(n):
		nsum += v[i]*w[i]
		ndenom += w[i]

	return nsum/ndenom



def mean(v):

	n = len(v)
	vsum = 0
	vmean = 0

	for i in range(n):
		vsum += v[i]

	vmean = vsum/float(n)

	return(vmean)

def geom_mean(v):

	n = len(v)
	lsum = 0

	if min(v) <= 0:
		return 0

	for i in range(n):
		lsum += log(v[i])

	lsum = lsum/n

	return exp(lsum)

def harmonic_mean(v):

	n = len(v)
	vsum = 0
	vmean = 0

	for i in range(n):
		vsum += 1.0/v[i]

	vmean = vsum/float(n)

	return(1.0/vmean)

#
# given the sim scores between two transcripts we can summarize it into a 
# single score with this function by producing a weighted mean based on the 
# lengths of the transcripts. 
# parameters:
# s1 is amount of t1 shared by t2
# s2 is amount of t2 shared by t1
def pair_sim(s1, s2, l1, l2):
	#
	# create weights from the lengths
	ll1 = log(l1)
	ll2 = log(l2)
	llmean = mean([ll1, ll2])
	w1 = 1/exp(ll1-llmean)
	w2 = 1/exp(ll2-llmean)

	rres = wmean([s1, s2], [w1, w2])

	return rres

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


#==============================================================================
# defs
#==============================================================================


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

def message(sz):
	sys.stderr.write("[seal-mappability-process] " + sz + "\n")


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Reduce the mappability table down to unique pairs")
parser.add_argument('file', type=str, help="Output of seal-mappability-graph.py")
parser.add_argument('gtf', type=str, help="GTF annotation for the transcriptome.")
parser.add_argument('fai', type=str, help="FAI index of the transcriptome containing feature lengths")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

