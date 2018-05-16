#!/usr/bin/python
#==============================================================================
# bam-gcounts-cs.py
#
# Shawn Driscoll
# 20170608
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Intersects bam alignments with an annotation. This is setup for 
# quantifying coordinate sorted alignments and handles quantification 
# in a per-locus fashion by reading through the alignments until a 
# gap is found, calling that a locus and sending the bundle of alignments 
# out to a worker function which then finds intersections with the annotation. 
#==============================================================================

import sys, argparse, math, re, os
from os.path import isfile, expanduser, isdir
from collections import defaultdict
from time import localtime
# used for printing exceptions
#import linecache
import traceback
import subprocess as sp
import hashlib
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock

# from igraph import *
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
HBIN = 16000

# bundle stuff
MIN_BUNDLE = 0
MAX_BUNDLE_GAP = 50

TMP_DIR = "_bam_gcounts"

# sam flags
SAM_PAIRED = 0x1
SAM_PROPERLY_PAIRED = 0x2
SAM_UNAL = 0x4
SAM_MATE_UNAL = 0x8
SAM_REVERSED = 0x10
SAM_MATE_REVERSED = 0x20
SAM_FIRST_MATE = 0x40
SAM_SECOND_MATE = 0x80
SAM_SECONDARY = 0x100

# sam alignment fields
SAM_QNAME = 0
SAM_FLAG = 1
SAM_RNAME = 2
SAM_POS = 3
SAM_MAPQ = 4
SAM_CIGAR = 5
SAM_RNEXT = 6
SAM_PNEXT = 7
SAM_TLEN = 8
SAM_SEQ = 9
SAM_QUAL = 10
SAM_FIRST_ATTR = 11
# after replacement
SAM_RLEN = 9
SAM_ALNLEN = 10

# region fields
REGION_RNAME = 0
REGION_START = 1
REGION_END = 2
REGION_STRAND = 3
REGION_TAG = 4
REGION_INDEX = 5

# GTF fields
GTF_RNAME = 0
GTF_SOURCE = 1
GTF_FEATURE = 2
GTF_START = 3
GTF_END = 4
GTF_SCORE = 5
GTF_STRAND = 6
GTF_FRAME = 7
GTF_ATTRIBUTE = 8

#==============================================================================
# main
#==============================================================================

def main_mp(args):

	max_gap = MAX_BUNDLE_GAP
	aln_buff = []
	dhits = defaultdict(float)
	phits = 0
	khits = 0
	# initial values
	frag_len = {"mean":200, "stdev":80**2, "n":1}

	num_parsed = 0
	num_frags = 0
	num_aligned = 0
	num_assigned = 0
	num_assigned_multi = 0
	num_passq = 0
	num_bundles = 0
	num_buff = 0
	num_no_bundle = 0
	num_buffered = 0
	num_primary = 0
	num_secondary = 0

	fout_stub = ""	
	if args.o is None:
		# setup output file stub
		tmp = (args.bam).split("/")
		# last part is the file name
		tmp2 = tmp[-1].split(".")
		fout_stub = tmp2[0]
	else:
		fout_stub = args.o
	
	# multiprocessing stuff
	tasks = JoinableQueue()
	results = Queue()
	p = None
	pool = []
	tlock = Lock()
	
	# make a stub name for all of the bundle outputs to disk
	tmp_name = hashlib.md5(args.bam).hexdigest()
	
	#
	# announce what we're processing
	#
	message("Processing {}".format(args.bam))
	
	#
	# parse the gtf and build lookup table
	#
	message("Loading {}".format(args.gtf))
	dannot, dgid2tid, dgname2tid = gtf_load(args.gtf)
	message("building index")
	lktable = gtf_lookup_table(dannot)
	message("building region index")
	dregions = bundle_chrom_regions(dannot)
	message("done")
	
	#
	# create temp folder for dumping parsed sam lines
	#
	if not isdir(TMP_DIR):
		# create temporary directory
		os.mkdir(TMP_DIR)
	
	r = None
	r0 = None
	bundle = ""
	bundleLast = ""
	
	#
	# start processes
	#
	for i in range(args.p):
		p = Process(target=bundle_worker, args=(tasks, results, args, dannot, tlock, ))
		p.daemon = True
		p.start()
		pool.append(p)
	
	message("processing alignments")
	fin = open_bam(args.bam)
	
	# open initial bundle output
	fout_bundle = "{}/{}.{}.sam".format(TMP_DIR, tmp_name, num_bundles)
	fout = open(fout_bundle, "w")
	
	for szl in fin:
		# skip header lines
		if szl[0]=="@":
			continue
		
		# parse
		num_parsed += 1
		aln = sam_parse(szl)
				
		if sam_unaligned(aln):
			continue
		
		num_aligned += 1
		
		if aln[SAM_MAPQ] < args.min_mapq:
			continue
		
		num_passq += 1
		
		sam_aln_extend(aln)
		
		start = sam_soft_start(aln)
		end = sam_end_pos(aln)
		r = region_init(aln[SAM_RNAME], start, end, sam_strand(aln))
	
		# check bundle
		rhits = bundle_lookup_hits(r, dregions)
		if len(rhits) != 1:
			num_no_bundle += 1
			continue
		
		bundle = rhits[0]
		
		# 
		# we could have a new bundle id OR we could have a null bundle id but a buffer with hits 
		# in it
		if bundle != bundleLast:
			# finished collecting hits from a bundle, time to process
			# the buffer...maybe...
			if bundleLast != "":
				progress_message("processing {}".format(bundleLast))

			if (not bundle_same_chrom(bundle, bundleLast)) or num_buffered > MIN_BUNDLE:
				# close buffer file and put into task queue
				fout.close()
				tasks.put(fout_bundle)
				# start a new one				
				num_bundles += 1
				fout_bundle = "{}/{}.{}.sam".format(TMP_DIR, tmp_name, num_bundles)
				fout = open(fout_bundle, "w")
			
				# clear alignment buffer count
				num_buffered = 0
		
		#
		# write current read to bundle file
		fout.write(szl)
		num_buffered += 1
		# keep bundle name
		bundleLast = rhits[0]
		
	fin.close()
	
	if num_buffered > 0:
		progress_message("processing {}".format(bundleLast))
		# close 
		fout.close()
		tasks.put(fout_bundle)

	# get out of the progress message line
	sys.stderr.write("\n")
	
	message("waiting for child processes to finish")
	
	for p in pool:
		tasks.put(None)
	
	tasks.join()
	for p in pool:
		p.join()
	
	results.put(None)
	
	# merge the individual buffer outputs
	message("merging child process outputs")
	
	dhits = defaultdict(float)
	frag_means = []
	frag_sds = []
	frag_ns = []
	while True:
		fname = results.get()
		if fname is None:
			break
		
		with open(fname, "r") as fin:
			szheader = fin.readline()
			# parse the header line
			aln = szheader.strip().split("\t")
			dheader = {}
			for k in aln:
				tmp = k.split(":")
				dheader[tmp[0]] = float(tmp[1])
			
			num_frags += dheader['num_frags']
			num_assigned += dheader['num_assigned']
			num_assigned_multi += dheader['num_assigned_multi']
			frag_means.append(float(dheader['frag_mean']))
			frag_sds.append(float(dheader['frag_stdev']))
			frag_ns.append(float(dheader['frag_n']))
			
			for szl in fin:
				aln = szl.strip().split("\t")
				dhits[aln[0]] += float(aln[1])
		
		# drop the process hits file
		os.unlink(fname)
	
	for i in range(len(frag_sds)):
		frag_sds[i] = frag_sds[i]/(frag_ns[i]-1)
	
	
	# finalze the fragment length calculation
	frag_len_mean = weighted_mean(frag_means, frag_ns)
	frag_len_stdev = math.sqrt(weighted_mean(frag_sds, frag_ns))	
	
	message("writing results to stdout")
	
	if args.exon:
		generate_exon_output(dannot, dhits, fout_stub)
	else:
		generate_transcript_output(dannot, dhits, frag_len_mean, fout_stub)
	
	message("#mean fragment: {:0.2f} +/- {:0.2f}".format(frag_len_mean, frag_len_stdev), False)
	message("#num parsed: {}".format(num_parsed), False)
	message("#num aligned: {}".format(num_aligned), False)
	message("#num passq: {}".format(num_passq), False)
	message("#num no bundle: {}".format(num_no_bundle), False)
	message("#num frags: {}".format(num_frags), False)
	message("#num assigned: {}".format(num_assigned), False)
	message("#num multi: {}".format(num_assigned_multi), False)
	message("#num buffered: {}".format(num_buff), False)
	message("#total processed bundles: {}".format(num_bundles), False)
	
	# remove the temp dir
	try:
		os.rmdir(TMP_DIR)
	except:
		warning_message("Could not remove temp folder {}. Maybe it isn't empty?".format(TMP_DIR))
		pass
	
	return 0

#
# bundle_worker
# worker function for child processes. this worker receives file names of 
# SAM data as tasks. for each task the file is read in and parsed 
# via 'sam_parse' and 'sam_aln_extend'. the alignments are bundled into 
# a list and passed off to 'process_bundle'. each instance of the worker
# maintains it's own hits table and lookup table to avoid collisions between
# processes.
def bundle_worker(task_queue, results_queue, args, annot, tlock):
	name = current_process().name
	
	# each worker needs its own lookup table
	lktable = gtf_lookup_table(annot)
	dhits = defaultdict(float)
	
	# initial values
	frag_len = {"mean":200, "stdev":80**2, "n":1}
	num_parsed = 0
	num_frags = 0
	num_aligned = 0
	num_assigned = 0
	num_frags = 0
	num_assigned_multi = 0
	num_passq = 0
	num_bundles = 0
	num_buff = 0
	num_no_bundle = 0	
	
	aln_buff = []
	
	while True:
		item = task_queue.get()
		if item is None:
			task_queue.task_done()
			break
		
		#message_mp("received {}".format(item), name, tlock)
		
		# read file
		aln_buff = []
		with open(item, "r") as fin:
			for szl in fin:
				aln = sam_parse(szl)
				sam_aln_extend(aln)
				aln_buff.append(aln)
		
		if len(aln_buff) > 0:
			# process
			tmp0, tmp1, tmp2 = process_buffer(aln_buff, args, lktable, annot, dhits, frag_len)
			num_frags += tmp0
			num_assigned += tmp1
			num_assigned_multi += tmp2

		#message_mp("finished {}".format(item), name, tlock)
		os.unlink(item)
		# finished task
		task_queue.task_done()
		
	
	# at the end we have to write this worker's hit results out to a file and 
	# send that file name back to 'main' via the 'results_queue'
	
	fout_name = "{}.hits".format(name)
	with open(fout_name, "w") as fout:

		# write a header line that passes back the fragment length info as well as 
		# other stats
		fout.write("frag_mean:{}\t".format(frag_len['mean']))
		fout.write("frag_stdev:{}\t".format(frag_len['stdev']))
		fout.write("frag_n:{}\t".format(frag_len['n']))
		fout.write("num_frags:{}\t".format(num_frags))
		fout.write("num_assigned:{}\t".format(num_assigned))
		fout.write("num_assigned_multi:{}\n".format(num_assigned_multi))
		
		for k in dhits.keys():
			fout.write("{}\t{}\n".format(k, dhits[k]))
	
	message_mp("wrote {}".format(fout_name), name, tlock)
	results_queue.put(fout_name)
	
	return

def main(args):

	fin = open_bam(args.bam)
	max_gap = 100
	aln_buff = []
	dhits = defaultdict(float)
	phits = 0
	khits = 0
	# initial values
	frag_len = {"mean":200, "stdev":80**2, "n":1}
	num_parsed = 0
	num_frags = 0
	num_aligned = 0
	num_assigned = 0
	num_assigned_multi = 0
	num_passq = 0
	num_bundles = 0
	num_buff = 0
	num_no_bundle = 0
	
	#
	# parse the gtf and build lookup table
	#
	message("Loading {}".format(args.gtf))
	dannot, dgid2tid, dgname2tid = gtf_load(args.gtf)
	message("building index")
	lktable = gtf_lookup_table(dannot)
	message("building region index")
	dregions = bundle_chrom_regions(dannot)
	message("done")
		
	r = None
	r0 = None
	bundle = ""
	bundleLast = ""
	
	message("processing alignments")
	for szl in fin:
		# skip header lines
		if szl[0]=="@":
			continue
		
		# parse
		num_parsed += 1
		aln = sam_parse(szl)
				
		if sam_unaligned(aln):
			continue
		
		num_aligned += 1
		
		if aln[SAM_MAPQ] < args.min_mapq:
			continue
		
		num_passq += 1
		
		sam_aln_extend(aln)
		
		start = sam_soft_start(aln)
		end = sam_end_pos(aln)
		r = region_init(aln[SAM_RNAME], start, end, sam_strand(aln))
	
		# check bundle
		rhits = bundle_lookup_hits(r, dregions)
		if len(rhits) != 1:
			num_no_bundle += 1
			continue
		
		bundle = rhits[0]
		
		# 
		# we could have a new bundle id OR we could have a null bundle id but a buffer with hits 
		# in it
		if bundle != bundleLast:
			# finished collecting hits from a bundle, time to process
			# the buffer
			if len(aln_buff) > 0:
				progress_message("processing {}".format(bundleLast))
				num_bundles += 1
				num_buff += len(aln_buff)
				
				tmp1, tmp2 = process_buffer(aln_buff, args, lktable, dannot, dhits, frag_len)
				num_assigned += tmp1
				num_assigned_multi += tmp2
			
			# clear alignment buffer
			aln_buff = []
		
		#
		# add current read to buffer
		aln_buff.append(aln)
		# keep bundle name
		bundleLast = rhits[0]
		
	fin.close()
	
	if len(aln_buff) > 0:
		progress_message("processing {}".format(bundleLast))
		num_buff += len(aln_buff)
		tmp1, tmp2 = process_buffer(aln_buff, args, lktable, dannot, dhits, frag_len)
		num_assigned += tmp1
		num_assigned_multi += tmp2
	
	total_hits = 0
	
	# get out of the progress message line
	sys.stderr.write("\n")
	
	if args.exon:
		generate_exon_output(dannot, dhits)
	else:
		generate_transcript_output(dannot, dhits, frag_len)
	
	message("#mean fragment: {:0.2f} +/- {:0.2f}".format(frag_len['mean'], math.sqrt(frag_len['stdev']*1.0/(frag_len['n']-1))))
	message("#num parsed: {}".format(num_parsed))
	message("#num aligned: {}".format(num_aligned))
	message("#num passq: {}".format(num_passq))
	message("#num no bundle: {}".format(num_no_bundle))
	message("#num assigned: {}".format(num_assigned))
	message("#total counts: {}".format(total_hits))
	message("#num multi: {}".format(num_assigned_multi))
	message("#num buffered: {}".format(num_buff))
	return 0

#==============================================================================
# general functions
#==============================================================================

def generate_exon_output(annot, dhits, fstub):
	
	# write results to file
	fout_name = "{}.exon_quant".format(fstub)
	message("Writing results to {}".format(fout_name))
	with open(fout_name, "w") as fout:
	
		head_out = "transcript_id\tgene_id\tgene_name\tchrom\tstrand\tnum_exons\tlengths\thits\taverage"
		fout.write(head_out + "\n")
		
		n = 0
		
		for tid in sorted(annot.keys()):
			# loop through exons and collect info
			lengths = []
			hits = []
			mean_hits = []
			n = len(annot[tid]['exons'])
			for i in range(n):
				elen = annot[tid]['exons'][i][1]-annot[tid]['exons'][i][0]+1
				lengths.append(elen)
				k = "{}:{}".format(tid, i)
				if k in dhits:
					hits.append(dhits[k])
				else:
					hits.append(0)
			
			# calculate average counts/depth at each exon
			for i in range(len(hits)):
				if lengths[i] > 0:
					mean_hits.append(hits[i]*1.0/lengths[i])
				else:
					mean_hits.append(0)
	
			lout = [
				tid, annot[tid]['gene_id'], annot[tid]['gene_name'], annot[tid]['rname'], annot[tid]['strand'], 
				"{}".format(len(hits)), 
				",".join(map(str, lengths)), 
				",".join(["{:0.4f}".format(x) for x in hits]), 
				",".join(["{:0.4f}".format(x) for x in mean_hits])
			]
			
			fout.write("\t".join(lout))
			fout.write("\n")
	
	return 0

def generate_transcript_output(annot, dhits, frag_len_mean, fstub):

	# write results to file
	fout_name = "{}.ghits".format(fstub)
	message("Writing results to {}".format(fout_name))

	head_out = "transcript_id\tgene_id\tgene_name\tchrom\tstrand\tlength\teff_length\thits\teff_hits\tfpkm\ttpm"
	dheader = { 'length':5, 'eff_length':6, 'hits':7, 'eff_hits':8, 'fpkm':9 }

	total_hits = 0
	total_fpkm = 0
	outlines = []
	# this first loop starts to build each output line and counts up total hits
	for tid in sorted(annot.keys()):
		#print "{}\t{}".format(tid, dhits[tid])
		lout = [tid]
		hits = 0
		if tid in dhits:
			hits = dhits[tid]
		eff_hits = hits

		length = annot[tid]['length']
		eff_length = annot[tid]['length']-frag_len_mean
		if eff_length <= 4:
			eff_length = 0
			eff_hits = 0
		else:
			eff_hits = float(hits)*annot[tid]['length']/eff_length
		
		lout += [annot[tid]['gene_id'], annot[tid]['gene_name'], annot[tid]['rname'], annot[tid]['strand']]
		
		lout += [length, eff_length, hits, eff_hits]
		total_hits += hits
		outlines.append(lout)
		
		#print "\t".join(map(str, lout))
	
	# make fpkm
	
	# NOTE: effective length is used to describe the effective transcript length
	# of the observed hit count. effective hits is used to describe the estimated
	# hits based on the full transcript length. so when calculating FPKM we should 
	# use either observed count and effective length or effective count and full 
	# length. 
	
	nlines = len(outlines)
	for i in range(nlines):
		f = 0
		ll = outlines[i]
		if ll[dheader['hits']] > 0 and ll[dheader['eff_length']] > 0:
			f = float(ll[dheader['hits']])*1e9/(ll[dheader['eff_length']]*total_hits)
		
		total_fpkm += f
		outlines[i].append(f)

	# calculate TPM and write results
	with open(fout_name, "w") as fout:
	
		fout.write(head_out)
		fout.write("\n")
		
		for i in range(nlines):
			ll = outlines[i]
			tpm = 0
			if ll[-1] > 0:
				tpm = ll[dheader['fpkm']]*1e6/total_fpkm
			
			ll.append(tpm)
			for i in range(dheader['eff_length'], len(ll)):
				# format the floats to 4 decimals
				tmp = "{:0.4f}".format(ll[i])
				ll[i] = tmp
			fout.write("\t".join(map(str, ll)))
			fout.write("\n")
	
	return 0

def bundle_rname(sz):
	tmp = sz.split(":")
	return tmp[0]

def bundle_same_chrom(b1, b2):
	rname1 = bundle_rname(b1)
	return rname1==bundle_rname(b2)

def bundle_lookup_hits(r, lktable):
	
	hits = set()
	
	# hash the region
	rhk = region_hash(r)
	for k in rhk:
		if k in lktable:
			for r0 in lktable[k]:
				rres = compare_regions(r, r0)
				if rres > 0:
					# hit it
					hits.add(r0[REGION_TAG])
	
	return list(hits)

#
# this function bundles annotated transcript boundaries by chrosome
# into a dict. and build a lookup table.
def bundle_chrom_regions(dannot):
	
	tid = None
	rname = None
	dbundles0 = defaultdict(list)
	dbundles = defaultdict(list)
	dRegionLktable = defaultdict(list)

	#
	# build dict by chromosome name. each chrom name will have a list of 
	# regions. the regions start out as just the start/end boundaries of 
	# each transcript. 
	for tid in dannot.keys():
		r = region_init(dannot[tid]['rname'], dannot[tid]['start'], dannot[tid]['end'], 
			dannot[tid]['strand'])
		
		dbundles0[dannot[tid]['rname']].append(r)
	
	#
	# sort and collapse overlapping transcripts into bundles
	for rname in dbundles0.keys():
		# check for more than one region. if there's just one then we 
		# can just add it in if not then we have to sort and collapse
		if len(dbundles0[rname]) > 1:
			# sort
			elsort = sorted(dbundles0[rname], key=region_sort_key)
			# now collapse
			r0 = list(elsort[0])
			for i in range(len(elsort)):
				r = elsort[i]
				# check overlap
				rres = compare_regions(r, r0)
				# also join near neighbors (< N bases)
				if rres > 0 or (r[REGION_START]-r0[REGION_END] < MAX_BUNDLE_GAP):
					# combine
					r0[REGION_START] = min([r[REGION_START], r0[REGION_START]])
					r0[REGION_END] = max([r[REGION_END], r0[REGION_END]])
				else:
					# no overlap, check gap
					dbundles[rname].append(list(r0))
					r0 = list(elsort[i])
					
			# append the last one
			dbundles[rname].append(list(r0))
			
		else:
			# only one region so just add it to the final set
			dbundles[rname].append(dbundles0[rname][0])
	
	for rname in dbundles.keys():
		for r in dbundles[rname]:
			hk = region_hash(r)
			rstr = region_str(r)
			r[REGION_TAG] = rstr
			for hkey in hk:
				dRegionLktable[hkey].append(list(r))
	
	return dRegionLktable
	

def process_buffer(lbuff, args, lktable, annot, dhits, frag_len):
	
	n = len(lbuff)
	qname = ""
	qname_last = ""
	qname_buff = []
	flen0 = 0
	fstd0 = 0
	phits = 0
	khits = 0
	
	num_frags = 0
	num_assigned = 0
	num_assigned_multi = 0

	if n > 1:
		#message("sorting alignments by read name")
		lbuff_sorted = sort_alignments(lbuff)
	else:
		lbuff_sorted = lbuff
	
	# loop through the alignments and bundle by read name
	for aln in lbuff_sorted:
		qname = aln[SAM_QNAME]
		
		# check alignment for hits
		rres = process_alignment(args, lktable, aln)
		foo = sam_hits(aln)
		
		if qname != qname_last:
			num_frags += 1
				
			if len(qname_buff) > 0:
				# deal with it
				rres, flen = process_read(qname_buff, annot, dhits, args, frag_len)

				if rres==0 and len(foo) > 0:
					khits += 1
					
				if rres > 0:
					# we hit something. update mean and stdev of fragment length
					phits += 1
					num_assigned += 1
					frag_len['n'] += 1
					flen0 = frag_len['mean']
					fstd0 = frag_len['stdev']
					frag_len['mean'] = flen0 + (flen-flen0)*1.0/frag_len['n']
					frag_len['stdev'] = fstd0 + (flen-flen0)*(flen-frag_len['mean'])
				
				if rres > 1:
					num_assigned_multi += 1
				
				qname_buff = []
		
		qname_last = qname
		qname_buff.append(aln)

	if len(qname_buff) > 0:
		# deal with it
		num_frags += 1
		rres, flen = process_read(qname_buff, annot, dhits, args, frag_len)
		
		if rres > 0:
			# we hit something. update mean and stdev of fragment length
			num_assigned += 1
			frag_len['n'] += 1
			flen0 = frag_len['mean']
			fstd0 = frag_len['stdev']
			frag_len['mean'] = flen0 + (flen-flen0)*1.0/frag_len['n']
			frag_len['stdev'] = fstd0 + (flen-flen0)*(flen-frag_len['mean'])

		if rres > 1:
			num_assigned_multi += 1
	
	return num_frags, num_assigned, num_assigned_multi


#
# check single alignment, aln, for valid hits to the GTF annotation. the 
# alignment's hit list is modified and nothing is returned from this function
def process_alignment(args, lktable, aln):

	# list for hits to this alignment
#	aln['hits'] = []
	hit_index = len(aln)-1
	tmp = []
	dfinal = {}
	# get alignment strand
	aln_strand = sam_strand(aln)
	# set min overlap ratio for full alignment
	min_ratio = args.min_overlap_ratio
	# get alignment regions
	aln_r = sam_parse_regions(aln)
	target_sense = None
	
	if args.fr_stranded:
		if sam_paired(aln) and sam_first_mate(aln):
			target_sense = False
		elif sam_paired(aln) and (not sam_first_mate(aln)):
			target_sense = True
		elif not sam_paired(aln):
			target_sense = False
	elif args.rf_stranded:
		if sam_paired(aln) and sam_first_mate(aln):
			target_sense = True
		elif sam_paired(aln) and not sam_first_mate(aln):
			target_sense = False
		elif not sam_paired(aln):
			target_sense = True
				
	#
	# check for hits in each region of the alignment
	for r in aln_r:
		
		has_hits, d = gtf_find_hits(lktable, r, target_sense)
		if has_hits:
			tmp.append(d)
		
	if args.exon:
		# exon level quantification does not require much work. we keep hits as they are
		if len(tmp) > 0:
			dfinal = tmp[0]
			if len(tmp) > 1:
				for i in range(1, len(tmp)):
					d = tmp[i]
					for tid in d.keys():
						# see if we need to combine these or what
						if tid not in dfinal:
							dfinal[tid] = d[tid]
						else:
							# collision. combine them into a single hit. this can happen when a 
							# read is spliced but still both ends of the splice hit the same
							# exon. there's probably also a case when the two ends hit separate
							# exons.
							dfinal[tid]['length'] += d[tid]['length']							
		
		# check exon level overlap. if we're doing depth then everything can pass through
		for tid in dfinal.keys():
			if args.depth:
				# good, assign it
				aln[hit_index].append([tid, dfinal[tid]['index'], dfinal[tid]['length']])
			else:
				# check --eovl
				if dfinal[tid]['length'] >= args.eovl:
					# good, assign it
					aln[hit_index].append([tid, dfinal[tid]['index'], dfinal[tid]['length']])

	else:

		# transcript level. each hit 'tid' has to have the exon index removed
		# from it to build a dict by transcript id
		dfinal = {}
		if len(tmp) > 0:
			for i in range(0, len(tmp)):
				d = tmp[i]
				for tid0 in d.keys():
					# see if we need to combine these or what
					tid = tag_tid(tid0)
					if tid not in dfinal:
						dfinal[tid] = d[tid0]
					else:
						# collision. combine them into a single hit
						dfinal[tid]['index'].update(d[tid0]['index'])
						dfinal[tid]['length'] += d[tid0]['length']

		#
		# evaluate hit overlap length relative to acceptable ratio
		for tid in dfinal.keys():
			if args.depth:
				# no need to evaluate overlap ratio
				aln[hit_index].append([tid, dfinal[tid]['index'], dfinal[tid]['length']])
			else:
				rat = dfinal[tid]['length']*1.0/aln[SAM_ALNLEN]
				if rat >= min_ratio:
					# append hit. transcript id and the exons that were hit
					aln[hit_index].append([tid, dfinal[tid]['index'], dfinal[tid]['length']])

	return 0

#
# rbuffer is a set of alignments all with the same name. this function 
# checks hits for the alignments, establishes the alignment weight and 
# assigns the hits to the elements in the passed results dict
def process_read(rbuffer, annot, hit_table, args, frag_len):

	pe = False
	result = 0
	hits = set()
	gene_hits = set()

	aln1 = [] 
	aln2 = []
	#frag_len = []
	use_frag_len = True
	dtmp = defaultdict(int)
	this_flen = []

	# sort the reads into left and right
	for aln in rbuffer:
		if sam_paired(aln):
			pe = True
			if sam_first_mate(aln):
				aln1.append(aln)
			else:
				aln2.append(aln)
		else:
			aln1.append(aln)

	if pe:
		# pair the alignments
		pares = []
		dtmp = {}
	
		for i in range(len(aln1)):
			aln = aln1[i]
			aln_se = [sam_soft_start(aln), sam_end_pos(aln)]
			
			for j in range(len(aln2)):
				mate = aln2[j]
				
				if mate is None:
					continue

				mate_se = [sam_soft_start(mate), sam_end_pos(mate)]
				
				if not sam_are_mates(aln, mate):
					continue
				
				# they are mates. since the alignment hits are in terms of 'transcript_id:exon_index'
				# to pair up the hits I have to group them by transcript. this is also because
				# in some cases a single alignment may hit the same transcript more than one time. 
				# in the 'bam-gcounts-dev3.py' version this isn't an issue because the hits are 
				# reduced to transcript id already in 'process_alignment'
				
				da1hits = defaultdict(list)
				for hit in sam_hits(aln):
					tid = tag_tid(hit[0])
					da1hits[tid].append(hit)

				da2hits = defaultdict(list)
				for hit in sam_hits(mate):
					tid = tag_tid(hit[0])
					da2hits[tid].append(hit)
				
				# pair up the hits. that is to say we only want to move forward with cases when
				# this mate-pair are both hitting the same transcript
				for tid1 in da1hits.keys():
					for tid2 in da2hits.keys():
						if tid1 != tid2:
							continue
						
						# match, we can sort out the insert size now and potentially
						# assign the hits

						#
						# figure out within transcript insert size. the z-score
						# of the insert size relative to the running fragment length 
						# mean and stdev will help guide the assignment of the read
						# to the transcript or not
						#
						
						start = min(aln_se+mate_se)
						end = max(aln_se+mate_se)
						
						if aln_se[0] < mate_se[0]:
							left = aln
							right = mate
						else:
							left = mate
							right = aln
						
						left_bounds = sam_parse_regions(left)[0]
						right_bounds = sam_parse_regions(right)[-1]
						
						#
						# find start/end coordinates in transcript coordinates
						tid = tid1
						estart = 0
						eend = 0
						n = len(annot[tid]['exons'])
						
						# 
						# figure out the start and end exons of the pair
						#
						
						found_start = False
						found_end = False
						
						while estart < n:
							e = annot[tid]['exons'][estart]
							if left_bounds[REGION_START] <= e[1] and left_bounds[REGION_END] >= e[0]:
								found_start = True
								break
							estart += 1
							
						eend = estart
						
						while eend < n:
							e = annot[tid]['exons'][eend]
							if right_bounds[REGION_START] <= e[1] and right_bounds[REGION_END] >= e[0]:
								found_end = True
								break
							eend += 1
						
						if found_start and found_end:
							# found both start and end so we're good. figure out alignment 
							# insert size in terms of the transcript
						
							if eend == estart:
								isize = end-start+1
							else:
								# translate it
								tstart = max([0, start - annot[tid]['exons'][estart][0]]) + annot[tid]['texons'][estart][0]
								tend = max([0, end - annot[tid]['exons'][eend][0]]) + annot[tid]['texons'][eend][0]
								isize = tend-tstart+1
							
							# use insert size to filter out inappropriate (emperically) insert sizes.
							# this helps cut down on mis assignment of reads to transcripts
							frag_stdev = math.sqrt(frag_len['stdev']/max([1, frag_len['n']-1]))									
							isize_z = (isize-frag_len['mean'])*1.0/frag_stdev
							if isize > 0 and isize_z < 5:
								# passes the insert size filter so we can assign the hits
								this_flen.append(isize)
								# assign to the transcript's exons
								if args.depth:
									# add depth to each exon
									for hit in da1hits[tid1]:
										dtmp[hit[0]] += hit[2]
										
									for hit in da2hits[tid2]:
										dtmp[hit[0]] += hit[2]
										
								else:
									
									# add a single hit per exon unless the exon has already by 
									# tagged by this fragment
									for hit in da1hits[tid1]:
										if hit[0] not in dtmp:
											dtmp[hit[0]] = 1
																		
									for hit in da2hits[tid1]:
										if hit[0] not in dtmp:
											dtmp[hit[0]] = 1
						
			
				# finished with processing of aln and 'mate'. we can 'None' the mate
				# position in aln2's list so that it is skipped next time
				aln2[j] = None

	else:
		#----------------------------------------------------------------------
		# 
		# single-end condition
		#
		#----------------------------------------------------------------------
		dtmp = {}
		for aln in aln1:
			# use aligned length as fragment length
			this_flen.append(aln[SAM_ALNLEN])
			# make list of the hits
			aln_hits = sam_hits(aln)
			if len(aln_hits) > 0:
				for hit in aln_hits:
					if hit[0] not in dtmp:
						if args.depth:
							dtmp[hit[0]] = hit[2]
						else:
							dtmp[hit[0]] = 1
					else:
						# this read hit the same transcript more than once.  that is to say the 
						# read is multi-mapped but within the same transcript. this is separate
						# an ALIGNMENT hitting the same transcript more than once which happens
						# all the time when an alignment is spliced. that is handled in 
						# 'process_alignment'
						if args.depth:
							# if we're counting depth take the maximum
							dtmp[hit[0]] = max([dtmp[hit[0]], hit[2]])
							# if we're counting hits then we only need to count this 
							# read at this feature 1 time. 
							warning_message("PR same read aligned to the same feature more than 1 time")
							
	# 
	# we can bail out here if there were no hits
	if len(dtmp.keys())==0:
		return 0, -1
	
	# make two lists: 1 is the list of targets and the other is a list of the values
	ltargets = []
	lweights = []
	result = 1
	for k in dtmp.keys():
		ltargets.append(k)
		lweights.append(dtmp[k])
	
	if args.exon:
		# at exon level we want to assign the full value to each element. 
		for i in range(len(ltargets)):
			hit_table[ltargets[i]] += lweights[i]
			
	else:
		# for transcript level we have to scale the weights as: w[i]/sum(w)*max(w)
		total_weight = sum(lweights)
		max_weight = max(lweights)
		
		num_genes = 1
		if args.G:
			# how many genes did we hit?
			for tid in ltargets:
				gene_hits.update([annot[tag_tid(tid)]['gene_id']])
			num_genes = len(list(gene_hits))
		
		if num_genes > 1:
			result = 2
		
		for i in range(len(ltargets)):
			wi = ((lweights[i]*1.0/total_weight)*max_weight) * 1.0/(num_genes**2)
			hit_table[ltargets[i]] += wi
	

	if len(this_flen) > 0:
		frag_mean = np.mean(this_flen)
	else:
		frag_mean = -1

	return result, frag_mean



def open_bam(fname):
	cmd = "samtools view {}".format(fname)
	p = sp.Popen(cmd.split(), stdout=sp.PIPE)
	return p.stdout


def weighted_mean(mu, wi):
	
	sum_prod = 0
	sum_wi = 0
	
	for i in range(len(mu)):
		sum_prod += mu[i]*wi[i]
		sum_wi += wi[i]
	
	return sum_prod/sum_wi

#
# for sorting
def aln_sort_name_key(item):
	return item[SAM_QNAME]

# 
# provided 'laln' is a list of parsed alignments			
def sort_alignments(laln):
	laln_sorted = sorted(laln, key=aln_sort_name_key)
	return laln_sorted

#==============================================================================
# GTF related functions
#==============================================================================


def gtf_parseline(sz):

	tmp = sz.strip().split("\t")

	grow = {
		'rname':tmp[0],
		'db':tmp[1],
		'type':tmp[2],
		'start':int(tmp[3]),
		'end':int(tmp[4]),
		'strand':tmp[6],
		'attrs':{}}

	# parse attributes
	fsplit = tmp[8].split("\"")
	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		grow['attrs'][key.strip()] = fsplit[i+1].strip()
		i += 2

	return grow

def gtf_transcript_id(grow):
	if "transcript_id" in grow['attrs']:
		return grow['attrs']['transcript_id']
	return None

def gtf_gene_id(grow):
	if "gene_id" in grow['attrs']:
		return grow['attrs']['gene_id']
	return None

def gtf_gene_name(grow):
	if "gene_name" in grow['attrs']:
		return grow['attrs']['gene_name']
	return None

def gtf_load(fname):
	#
	# load a gtf file into a few annotation tables

	dannot = {}
	dgid2tid = defaultdict(set)
	dgname2tid = defaultdict(set)

	fin = open(fname, "r")
	for szl in fin:
		grow = gtf_parseline(szl)

		tid = gtf_transcript_id(grow)
		if tid not in dannot:
			dannot[tid] = {
				'rname':grow['rname'],
				'strand':grow['strand'],
				'db':grow['db'],
				'gene_id':gtf_gene_id(grow),
				'gene_name':gtf_gene_name(grow),
				'num_exons':0,
				'start':grow['start'],
				'end':grow['end'],
				'length':0,
				'exons':[]}

		# update feature length
		dannot[tid]['length'] += grow['end']-grow['start']+1

		dannot[tid]['exons'].append([grow['start'], grow['end']])
		if grow['start'] < dannot[tid]['start']:
			dannot[tid]['start'] = grow['start']
		if grow['end'] > dannot[tid]['end']:
			dannot[tid]['end'] = grow['end']

		if dannot[tid]['gene_id'] is not None:
			dgid2tid[dannot[tid]['gene_id']].add(tid)

		if dannot[tid]['gene_name'] is not None:
			dgname2tid[dannot[tid]['gene_name']].add(tid)

	fin.close()

	#
	# pass through the annotation and sort the exons by position
	for tid in dannot.keys():
		dannot[tid]['exons'].sort(key=lambda x: x[0])
		# create transcript coordinate version of the exons
		dannot[tid]['texons'] = []
		left = 1
		right = left
		for e in dannot[tid]['exons']:
			l = e[1]-e[0]+1
			right = left + l - 1
			dannot[tid]['texons'].append([left, right])
			left = right+1
		
	return dannot, dgid2tid, dgname2tid

#
# build a quick-lookup table for exon features from a loaded gtf
def gtf_lookup_table(annot):

	lktable = defaultdict(list)

	# loop through annotation by transcript id
	for tid in annot.keys():
		# loop through exons
		eidx = 0
		for e in annot[tid]['exons']:
			# tags will be 'transcript_id:exon_index'
			r = region_init(annot[tid]['rname'], e[0], e[1], annot[tid]['strand'], tag="{}:{}".format(tid, eidx), index=eidx)
			eidx += 1
			h = region_hash(r)
			for hid in h:
				# insert the region in each binf
				lktable[hid].append(r)

	return lktable

def tag_tid(sz):
	tmp = sz.split(":")
	return tmp[0]

def tag_exonid(sz):
	tmp = sz.split(":")
	return int(tmp[1])

#
# look up region r in the lookup table. return info from each hit including 
# the tag and index values
def gtf_find_hits(lktable, r, target_sense):

	# hash the region, r
	lhash = region_hash(r)
	# hits list
	d_hits = {}
	has_hits = False

	# scan through hash
	for h in lhash:
		if h in lktable:
			# scan through regions in this bucket
			for r0 in lktable[h]:
				rres = compare_regions(r, r0)
				if rres > 0:
					sense = r[REGION_STRAND]==r0[REGION_STRAND]
					
					if (target_sense is None) or (target_sense == sense):

						# we have a hit!
						has_hits = True
						ovl_len = region_overlap_length(r, r0)
				
						if ovl_len[0] == 0:
							error_message("found overlap with length 0!")
							print region_str(r), region_str(r0)
							sys.exit(1)
	
						if r0[REGION_TAG] not in d_hits:
							try:
								d_hits[r0[REGION_TAG]] = { 'index': set([r0[REGION_INDEX]]), 'length': ovl_len[0] }
							except:
								print r0[REGION_TAG], r0[REGION_INDEX], ovl_len
								sys.exit(1)

					# I'm pretty sure the 'else' condition here is impossible. we are looking up hits
					# for a single alignment region. that single region could never hit more than 
					# one exon of any single transcript because no two exons of a single transcript
					# occupy the same space

#					else:
#						# add additional hit only if it is not to the same exon as one already
#						# recorded. this would happen when we have a region that crosses 
#						# the bin boundary and therefore would be found to hit the same feature
#						# if that feature also crosses the bin boundary
#						if r0['index'] not in set(d_hits[r0['tag']]['index']):
#							d_hits[r0['tag']]['index'].append(r0['index'])
#							d_hits[r0['tag']]['length'].append(ovl_len[0])

	#
	# return the dict of hits
	return has_hits, d_hits

#==============================================================================
# region related functions
#==============================================================================

def region_sort_key(r):
	return r[REGION_START]

def region_init(rname, start, end, strand, tag=None, index=None):
	# region is just a list
	r = [rname, int(start), int(end), strand, None, None]
	# set optional fields
	if tag is not None:
		r[REGION_TAG] = tag
	if index is not None:
		r[REGION_INDEX] = index
	
	return r

#
# assuming 'r' is a dict with 'rname', 'start' and 'end' fields
def region_hash(r):
	bin0 = binN = 0

	bin0 = int(r[REGION_START])/HBIN
	binN = int(r[REGION_END])/HBIN

	hout = ["{}:{}".format(r[REGION_RNAME], bin0)]
	
	if binN > bin0:
		while bin0 < binN:
			bin0 += 1
			hout.append("{}:{}".format(r[REGION_RNAME], bin0))

	return hout

#
# compare two regions. 
# return value:
# 0 for no overlap
# 1 for overlap
# 2 for identical
def compare_regions(r1, r2):

	rres = 0

	# check ref names. if not equal then we're done
	if r1[REGION_RNAME] != r2[REGION_RNAME]:
		return 0

	# ref names must be equal
	if r1[REGION_START]==r2[REGION_START] and r1[REGION_END]==r2[REGION_END]:
		# starts and ends are identical
		return 2

	# now check for overlap
	if r1[REGION_END] >= r2[REGION_START] and r2[REGION_END] >= r1[REGION_START]:
		# overlap!
		return 1
	
	return rres
	

def region_str(r):
	sz_pos = "{}:{}-{}".format(r[REGION_RNAME], r[REGION_START], r[REGION_END])

	if r[REGION_TAG] is not None:
		sz_pos += "|{}".format(r[REGION_TAG])

	if r[REGION_INDEX] is not None:
		sz_pos += "|{}".format(r[REGION_INDEX])

	return sz_pos

def region_length(r):
	return r[REGION_END]-r[REGION_START]+1

#
# calculate length of overlap between two regions.
# this is accomplished by finding the minimum value
# of 4 different measurements:
# A: length of r1
# B: length of r2
# C: end of r1 - start of r2
# D: end of r2 - start of r1
# if the minimum is negative then the result is 0
def region_overlap_length(r1, r2):
	len_A = region_length(r1)
	len_B = region_length(r2)
	len_C = r1[REGION_END]-r2[REGION_START]+1
	len_D = r2[REGION_END]-r1[REGION_START]+1

	rres = min([len_A, len_B, len_C, len_D])

	if rres <= 0:
		return [0, 0, 0]

	# return length of overlap as well as ratios of the overlap to the length 
	# of the features
	return [ rres, rres*1.0/len_A, rres*1.0/len_B ]



#==============================================================================
# SAM related functions
#==============================================================================

# if we have hits they will be in the last field of the alignment
def sam_hits(aln):
	return aln[-1]

#
# adjust start position for softclipping
def sam_soft_start(aln):
	cigar = aln[SAM_CIGAR]
	r = re.search("^([0-9]+)[S]", cigar)
	pos = aln[SAM_POS]
	if r:
		pos += int(r.group(1))
	return pos

#
# calculate end coordinate of the alignment
def sam_end_pos(aln):
	pos = sam_soft_start(aln)
	lens = re.findall("([0-9]+)[MDN]", aln[SAM_CIGAR])
	for l in lens:
		pos += int(l)
	
	return(pos-1)	

# expecting that passed object is the list that was split from a single
# sam alignment line
def sam_aligned_length(laln):
	cigar = laln[SAM_CIGAR]
	r = re.findall("([0-9]+)[MD]", cigar)
	return sum(map(int, r))

def sam_strand(laln):
	flag = int(laln[SAM_FLAG])
	s = "+"
	if flag & SAM_REVERSED:
		s = "-"
	return s

# return True if a pait of alignments are mates. this ignores the 0x40 and 0x80 
# flags because this program pre-sorts the alignments out into first/second 
# mate buckets. no need to double check
def sam_are_mates(laln1, laln2):
	
	# extra condition redundant for this code
	# if sam_first_mate(laln1) and not sam_first_mate(laln2):
	
	if laln1[SAM_RNEXT]==laln2[SAM_RNEXT]:
		if sam_strand(laln1) != sam_strand(laln2):
			if (int(laln1[SAM_PNEXT])==int(laln2[SAM_POS])) and (int(laln1[SAM_POS])==int(laln2[SAM_PNEXT])):
				return True
	
	return False 
	
def sam_header_parse(sz):
	ll = szl.strip().split("\t")
	d = {}
	d['type'] = re.sub("^@", "", ll[0])
	# parse additional fields
	for i in range(1, len(ll)):
		tmp = ll[i].split(":")
		d[tmp[0]] = tmp[1]
	
	return d

# 
# split a sam line
def sam_parse(sz):
	aln = sz.strip().split("\t")

	for i in [SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_MAPQ, SAM_TLEN]:
		aln[i] = int(aln[i])
	
	# replace the sesquence with the length of the sequence
	aln[SAM_SEQ] = len(aln[SAM_SEQ])
	# replace the quals with the aligned length
	aln[SAM_QUAL] = sam_aligned_length(aln)
		
	return aln

def sam_aln_extend(aln):
	# append a final slot that we can use to store hits
	aln.append(None) # this acts as a stop flag for functions that may loop through the alignment fields
	aln.append([]) # append hits list
	return 0	

def sam_parse_regions(aln):
	cigar = aln[SAM_CIGAR]
	lregions = []
	left = aln[SAM_POS]
	# operation index
	idx = 0
	# region index
	ridx = 0

	# extract all cigar operations. remember only M, D and N
	# actually advance anything. We also have to move forward if there
	# were any soft-clips
	op_len = [int(x) for x in re.findall("([0-9]+)[MIDNSHP\=X]", cigar)]
	op_type = re.findall("[0-9]+([MIDNSHP\=X])", cigar)
	
	if op_type[0]=="S":
		left += op_len[0]
		idx = 1

	# start right off at the left position
	right = left
	
	# loop through
	while idx < len(op_len):
		if op_type[idx] == "M" or op_type[idx] == "D":
			right += op_len[idx]
		elif op_type[idx] == "N":
			# found a junction. complete the previous region
			r = [aln[SAM_RNAME], left, right-1, sam_strand(aln), aln[SAM_QNAME], ridx]
			lregions.append(r)
			ridx += 1
			left = right + op_len[idx]
			right = left
			# move on...
		
		idx += 1
	
	# append region to list
	r = [aln[SAM_RNAME], left, right-1, sam_strand(aln), aln[SAM_QNAME], ridx]
	lregions.append(r)
	return lregions

#
# parse key/value paired attributes out from the parsed sam alignment
def sam_parse_atrributes(aln):
	# get the first position of the attributes
	idx = SAM_FIRST_ATTR
	dattr = {}
	
	while idx < len(aln):
		# loop until we hit the None entry that's added when we parse the
		# alignment into a list.
		if aln[idx] is None:
			break
		kv = aln[idx].split(":")
		dattr[kv[0]] = dattr[kv[2]]
		if kv[1]=="i":
			dattr[kv[0]] = int(dattr[kv[0]])
		elif kv[1]=="f":
			dattr[kv[0]] = float(dattr[kv[0]])
		
		idx += 1
	
	return dattr

#
# functions to check the status of the alignment flag
def sam_paired(aln):
	return ((aln[SAM_FLAG] & 0x1) != 0)

def sam_properly_paired(aln):
	return ((aln[SAM_FLAG] & 0x2) != 0)

def sam_unaligned(aln):
	return ((aln[SAM_FLAG] & 0x4) != 0)

def sam_mate_unaligned(aln):
	return ((aln[SAM_FLAG] & 0x8) != 0)

def sam_reversed(aln):
	return ((aln[SAM_FLAG] & 0x10) != 0)

def sam_mate_reversed(aln):
	return ((aln[SAM_FLAG] & 0x20) != 0)

def sam_first_mate(aln):
	return ((aln[SAM_FLAG] & 0x40) != 0)

def sam_secondary(aln):
	return ((aln[SAM_FLAG] & 0x100) != 0)

#
#
#

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))
	return 0

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))
	return 0

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))
	return 0

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0

def message_mp(sz, name, lock):
	lock.acquire()
	sys.stderr.write("[{}] {}\n".format(name, sz))
	lock.release()
	return 0

#def PrintException():
#    exc_type, exc_obj, tb = sys.exc_info()
#    f = tb.tb_frame
#    lineno = tb.tb_lineno
#    filename = f.f_code.co_filename
#    linecache.checkcache(filename)
#    line = linecache.getline(filename, lineno, f.f_globals)
#    sys.stderr.write('EXCEPTION IN ({}, LINE {} "{}"): {}\n'.format(filename, lineno, line.strip(), exc_obj))
#    return

def print_exception():
	exc_type, exc_value, exc_traceback = sys.exc_info()
	print "*** print_tb:"
	traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
	print "*** print_exception:"
	traceback.print_exception(exc_type, exc_value, exc_traceback,
	                          limit=2, file=sys.stdout)
	print "*** print_exc:"
	traceback.print_exc()
	print "*** format_exc, first and last line:"
	formatted_lines = traceback.format_exc().splitlines()
	print formatted_lines[0]
	print formatted_lines[-1]
	print "*** format_exception:"
	print repr(traceback.format_exception(exc_type, exc_value,
	                                      exc_traceback))
	print "*** extract_tb:"
	print repr(traceback.extract_tb(exc_traceback))
	print "*** format_tb:"
	print repr(traceback.format_tb(exc_traceback))
	print "*** tb_lineno:", exc_traceback.tb_lineno	
	return 0

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
#parser.add_argument('gtf', nargs="+", metavar="wig", type=str, 
#	help="Input wig or wig files to process")

parser.add_argument('gtf', type=str, 
	help="Annotation to quantify")
parser.add_argument('bam', type=str, 
	help="BAM alignments that are coordinate sorted.")

parser.add_argument('-p', action="store", default=1, type=int, 
	help="Number of processes to use for finding hits [1]")

parser.add_argument('-o', type=str, default=None, help="Stub for output file. Default is to use the input BAM file name.")

parser.add_argument('-G', action="store_const", const=True, default=False, 
	help="Enable multi-gene hit downweighting")

parser.add_argument('--exon', action="store_const", const=True, default=False, 
	help="Return exon-level quantification rather than transcript level")

parser.add_argument('--eovl', action="store", default=1, type=int, 
	help="Minimum number of bases of overlap for exon level counting [1]")

parser.add_argument('--depth', action="store_const", const=True, default=False, 
	help="Return depth per feature instead of read counts")

parser.add_argument("-l", "--min-overlap-ratio", default=0.96, type=float, action="store", 
	help="Minimum overlap as a ratio of the read, or read segment, length [0.96]")

parser.add_argument('-q', '--min-mapq', type=int, default=10, 
	help="Minimum MAPQ for counting alignment [10]")

#parser.add_argument('--depth-counts', action="store_const", const=True, default=False, 
#	help="Produce depth and average depth instead of read counts.")

#parser.add_argument('--exon-quant', action="store_const", const=True, default=False, 
#	help="Produce exon-level quantification instead of transcript level")

strand_group = parser.add_mutually_exclusive_group(required=False)
strand_group.add_argument('--fr-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a reverse stranded library [off]")

strand_group.add_argument('--rf-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a forward stranded library [off]")


args = parser.parse_args()

# 
# confirm input files exist
#
if not isfile(args.bam):
	error_message("input file {} does not exist".format(args.bam))
	sys.exit(1)

if not isfile(args.gtf):
	error_message("GTF file {} does not exist".format(args.gtf))
	sys.exit(1)

if __name__ == "__main__":

	try:
		sys.exit(main_mp(args))
#	except KeyboardInterrupt:
#		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)
