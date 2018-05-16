#!/usr/bin/python
#==============================================================================
# scrna-gcounts-bd.py
#
# Shawn Driscoll
# 20170614
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# With coordinate sorted alignments from a single well reads are 
# processed by annotated loci. We count MIs per gene and also 
# attempt to cluster MI's that are likely to have originated from the 
# same molecule but have diverged due to sequencing errors. 
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
#import igraph as ig

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

TMP_DIR = "_scrna_gcounts"

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

def main(args):

	aln_buff = []
	
	# count raw mi, corrected mi and reads
	dmi_hits = defaultdict(int)
	dmi_chits = defaultdict(int)
	dr_hits = defaultdict(int)
	
	phits = 0
	khits = 0
	# initial values
	num_parsed = 0
	num_secondary = 0
	num_frags = 0
	num_aligned = 0
	num_assigned = 0
	num_assigned_multi = 0
	num_passq = 0
	num_bundles = 0
	num_buff = 0
	num_sec_buff = 0
	num_no_bundle = 0
	num_multi_bundle = 0

	dmapq_buff = defaultdict(int)
	
	fout_stub = ""	
	if args.o is None:
		# setup output file stub
		tmp = (args.bam).split("/")
		# last part is the file name
		tmp2 = tmp[-1].split(".")
		fout_stub = tmp2[0]
	else:
		fout_stub = args.o

	message("BD Precise scRNA-seq quantification from alignments")	
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

	fin = open_bam(args.bam)
		
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
		
		if sam_secondary(aln):
			num_secondary += 1
		
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
			if len(rhits)==0:
				num_no_bundle += 1
			elif len(rhits) > 1:
				num_multi_bundle += 1
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
				#num_buff += len(aln_buff)
				
				tmp0, tmp1, tmp2 = process_buffer(aln_buff, args, lktable, dannot, dmi_hits, dmi_chits, dr_hits)
				num_frags += tmp0
				num_assigned += tmp1
				num_assigned_multi += tmp2
			
				# clear alignment buffer
				aln_buff = []
		
		#
		# add current read to buffer
		aln_buff.append(aln)
		num_buff += 1
		if sam_secondary(aln):
			num_sec_buff += 1
			
		dmapq_buff[str(aln[SAM_MAPQ])] += 1
		
		# keep bundle name
		bundleLast = rhits[0]
		
	fin.close()
	
	if len(aln_buff) > 0:
		progress_message("processing {}".format(bundleLast), last=True)
		#num_buff += len(aln_buff)
		tmp0, tmp1, tmp2 = process_buffer(aln_buff, args, lktable, dannot, dmi_hits, dmi_chits, dr_hits)
		num_frags += tmp0
		num_assigned += tmp1
		num_assigned_multi += tmp2
	
	message("Done. Writing results to {}.schits.".format(fout_stub))
	
	# write results to stdout
	generate_gene_output(dgid2tid, dannot, dmi_hits, dmi_chits, dr_hits, fout_stub)
	
	num_primary = num_aligned-num_secondary
	num_unal = num_parsed-num_aligned
	total_frags = num_primary+num_unal
	aligned_percent = num_primary*100.0/total_frags
	assigned_percent = num_assigned*100.0/num_frags
	
	message("# Results:", show_time=False)
	message("# parsed:               {}".format(num_parsed), show_time=False)
	message("# aligned:              {}".format(num_aligned), show_time=False)
	message("# secondary:            {}".format(num_secondary), show_time=False)
	message("# total frags:          {}".format(total_frags), show_time=False)
	message("# total frags aligned:  {} ({:0.2f}%)".format(num_primary, aligned_percent), show_time=False)
	message("# pass MAPQ:            {}".format(num_passq), show_time=False)
	message("# no bundle:            {}".format(num_no_bundle), show_time=False)
	message("# multi-bundle:         {}".format(num_multi_bundle), show_time=False)
	message("# processed:            {}".format(num_buff), show_time=False)
	message("# secondary processed:  {}".format(num_sec_buff), show_time=False)	
	message("# frags:                {}".format(num_frags), show_time=False)
	message("# assigned:             {} ({:0.2f}%)".format(num_assigned, assigned_percent), show_time=False)
	message("# assigned multi:       {}".format(num_assigned_multi), show_time=False)
	sys.stderr.write("\n")
	
	# generate log output
	with open("{}.runlog.tsv".format(fout_stub), "w") as fout:
		fout.write("# scrna-gcounts-bd.py\n")
		fout.write("# {}\n".format(time_string()))
		fout.write("# {}\n".format(args.bam))
		fout.write("#\n")
		fout.write("# Results:\n")
		fout.write("#stat\tcount\n")
		fout.write("num parsed\t{}\n".format(num_parsed))
		fout.write("num aligned\t{}\n".format(num_aligned))
		fout.write("num secondary\t{}\n".format(num_secondary))
		fout.write("total frags\t{}\n".format(total_frags))
		fout.write("total frags aligned\t{}\n".format(num_primary, aligned_percent))
		fout.write("num passq\t{}\n".format(num_passq))
		fout.write("num no bundle\t{}\n".format(num_no_bundle))
		fout.write("num multi-bundle\t{}\n".format(num_multi_bundle))
		fout.write("num processed\t{}\n".format(num_buff))
		fout.write("num secondary processed\t{}\n".format(num_sec_buff))
		fout.write("num frags\t{}\n".format(num_frags))
		fout.write("num assigned\t{}\n".format(num_assigned, assigned_percent))
		fout.write("num multi\t{}\n".format(num_assigned_multi))
		fout.write("\n")
		fout.write("# MAPQ histogram for processed reads\n")
		
		# create mapq report
		lmapq = sorted(map(int, dmapq_buff.keys()))
		fout.write("#MAPQ\tcount\n")
		for q in lmapq:
			k = str(q)
			fout.write("{}\t{}\n".format(k, dmapq_buff[k]))
		
	
	return 0

#==============================================================================
# MULTI-PROCESS main function (NOT CURRENTLY USED)
#==============================================================================

def main_mp(args):

	max_gap = MAX_BUNDLE_GAP
	aln_buff = []
	dhits = defaultdict(float)
	phits = 0
	khits = 0
	# initial values
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
		
		# skip unaligned or secondary hits
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
		generate_exon_output(dannot, dhits)
	else:
		generate_transcript_output(dannot, dhits, frag_len_mean)
	
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
			tmp0, tmp1, tmp2 = process_buffer(aln_buff, args, lktable, annot, dhits)
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
		fout.write("num_frags:{}\t".format(num_frags))
		fout.write("num_assigned:{}\t".format(num_assigned))
		fout.write("num_assigned_multi:{}\n".format(num_assigned_multi))
		
		for k in dhits.keys():
			fout.write("{}\t{}\n".format(k, dhits[k]))
	
	message_mp("wrote {}".format(fout_name), name, tlock)
	results_queue.put(fout_name)
	
	return


#==============================================================================
# general functions
#==============================================================================

def generate_gene_output(gid2tid, annot, dmi_hits, dmi_chits, dr_hits, fstub):
	
	# open output file
	with open("{}.schits".format(fstub), "w") as fout:
	
		head_out = "gene_id\tgene_name\tchrom\tstrand\tread_count\traw_mi_count\tmi_count"
		fout.write(head_out + "\n")
		
		# this first loop starts to build each output line and counts up total hits
		for gid in sorted(gid2tid.keys()):
			#print "{}\t{}".format(tid, dhits[tid])
			lout = [gid]
			tid = list(gid2tid[gid])[0]
			
			mi_hits = 0
			cmi_hits = 0
			r_hits = 0
	
			if gid in dmi_hits:
				mi_hits = dmi_hits[gid]
			if gid in dmi_chits:
				cmi_hits = dmi_chits[gid]
			if gid in dr_hits:
				r_hits = dr_hits[gid]
			
			lout += [annot[tid]['gene_name'], annot[tid]['rname'], annot[tid]['strand']]
			
			lout += [r_hits, mi_hits, cmi_hits]
			
			fout.write("\t".join(map(str, lout)) + "\n")
		
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
	
#
# this function takes a list of reads and checks them for hits to the
# annotation. if the list contains more than 1 read then the reads
# are sorted by read name
def process_buffer(lbuff, args, lktable, annot, dmi_hits, dmi_chits, dr_hits):
	
	n = len(lbuff)
	qname = ""
	qname_last = ""
	qname_buff = []
	phits = 0
	khits = 0
	
	num_frags = 0
	num_assigned = 0
	num_sec_assigned = 0
	num_assigned_multi = 0
	
	buff_hits = {}

	if n > 1:
		#message("sorting alignments by read name")
		lbuff_sorted = sort_alignments(lbuff)
	else:
		lbuff_sorted = lbuff
	
	# loop through the alignments and bundle by read name
	for aln in lbuff_sorted:
		qname = sam_qname(aln)
		
		# get hits for this alignment
		rres = process_alignment(args, lktable, aln)
		
		if (qname != qname_last) and len(qname_buff) > 0:
			
			num_frags += 1
				
			# deal with it
			rres, rhits = process_read(qname_buff)
			# rhits is a dict of gene ids and a single MI at each one
			
			if rres == 1:
				num_assigned += 1
				# this read hit a single gene. update the gene in 'buff_hits'
				# with this read's MI tags
					
				for k in rhits.keys():
						
					if k not in buff_hits:
						# add gene id dict. keep track of the count of each MI
						buff_hits[k] = {}
					
					# check for MI
					if rhits[k] not in buff_hits[k]:
						buff_hits[k][rhits[k]] = 0
					
					# increment count of this mi						
					buff_hits[k][rhits[k]] += 1
					# increment gene's read count
					dr_hits[k] += 1

			if rres > 1:
				num_assigned_multi += 1
			
			qname_buff = []
		
		qname_last = qname
		qname_buff.append(aln)

	# handle final assignment
	if len(qname_buff) > 0:
		num_frags += 1
		# deal with it
		rres, rhits = process_read(qname_buff)

		if rres == 1:
			num_assigned += 1
			# this read hit a single gene. update the gene in 'buff_hits'
			# with this read's MI tags
			for k in rhits.keys():

				if k not in buff_hits:
					# add gene id dict. keep track of the count of each MI
					buff_hits[k] = {}
				
				# check for MI
				if rhits[k] not in buff_hits[k]:
					buff_hits[k][rhits[k]] = 0
				
				# increment count of this mi						
				buff_hits[k][rhits[k]] += 1
				# increment gene's read count
				dr_hits[k] += 1


		if rres > 1:
			num_assigned_multi += 1

	# now we have to deal with the MI tags at each gene. since we won't be coming back to
	# this gene due to the nature of how this processing works we can finish it up now
	
	for gid in buff_hits.keys():
		mi_list = []
		mi_depth = []
		
		for mi in buff_hits[gid].keys():
			mi_list.append(mi)
			mi_depth.append(buff_hits[gid][mi])
		
		num_hits = len(mi_list)
		
		dmi_hits[gid] += num_hits
		
		if len(mi_list) > 1:
			#num_hits = mi_clustered_count(mi_list, mi_depth, thresh=1)
			num_hits = mi_cor_count(mi_list, mi_depth, thresh=1)
	
		dmi_chits[gid] += num_hits
	
	return num_frags, num_assigned, num_assigned_multi


#
# check single alignment, aln, for valid hits to the annotation. the 
# alignment's hit list is modified and nothing is returned from this function
def process_alignment(args, lktable, aln):

	# list for hits to this alignment
#	aln['hits'] = []
	hit_index = len(aln)-1
	tmp = []
	hits = set()
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

	# transcript level. each hit 'tid' has to have the exon index removed
	# from it to build a dict by transcript id
	dfinal = {}
	if len(tmp) > 0:
		for i in range(0, len(tmp)):
			d = tmp[i]
			for gid in d.keys():
				# see if we need to combine these or what
				if gid not in dfinal:
					dfinal[gid] = d[gid]
				else:
					# collision. combine them into a single hit
					dfinal[gid]['length'] += d[gid]['length']

	#
	# evaluate hit overlap length relative to acceptable ratio
	for gid in dfinal.keys():
		rat = dfinal[gid]['length']*1.0/aln[SAM_ALNLEN]
		if rat >= min_ratio:
			# append hit. transcript id and the exons that were hit
			hits.add(gid)
	
	#
	# each alignment must hit only one gene id or else it is ambiguous. due to the 
	# nature of the MI counting I have to throw these ones out. 
	hits = list(hits)
	if len(hits) > 1:
		# no hits
		aln[hit_index] = []
	else:
		aln[hit_index] = hits

	return 0

#
# aln1 is a list of alignments of the same read. here we parse out 
# the MI from the read name and assign it to the alignment's associated
# gene id in the 'hit_table'
def process_read(aln1):

	result = 0
	
	# keep a temporary dict of genes hit along with the MI tags from the read
	dtmp = {}

	for aln in aln1:
		# make list of the hits
		aln_hits = sam_hits(aln)
		# alignments come back with either 1 hit or none
		if len(aln_hits) == 1:
			gid = aln_hits[0]
			if gid not in dtmp:
				dtmp[gid] = sam_mi_from_qname(aln)
	
	# done for now
	nhit = len(dtmp.keys())
	if nhit > 0:
		result = 1
	if nhit > 1:
		result = 2

	return result, dtmp


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

#
# cluster MI tags using a graph with some thresholded maximum
# difference between mi strings to create edges
def mi_clustered_count(mi_list, mi_depth, thresh=1):
	
	n = len(mi_list)
	dmi = {}
	
	g = ig.Graph(directed=False)
	g.add_vertices(n)
	
	for i in range(n):
		dmi[mi_list[i]] = i
		
	edges = []
	edges_sz = []
	
	for i in range(n-1):
		for j in range(i+1, n):
			# compare
			diff = mi_hamming_dist(mi_list[i], mi_list[j])
			if diff <= thresh:
				# check depth. one of the two should appear to be a parent based
				# on depth
				depths = [mi_depth[i], mi_depth[j]]
				if max(depths) > 2*min(depths):
#					edges_sz.append((mi_list[i], mi_list[j]))
					edges.append((dmi[mi_list[i]], dmi[mi_list[j]]))
		
	if len(edges) > 0:
		g.add_edges(edges)
		# cluster
		gci = g.clusters(mode=ig.WEAK)
		n = len(gci)
		
	return n

#
# this more closely resembles the method used by BD. the idea is that 
# errors should be rare in the MIs so when combining, even with only
# a hamming dist of 1, there should be a parent with > 2*N counts.
def mi_cor_count(mi_list, mi_depth, thresh=1):
	# combine mi and depth
	mid = []
	for i in range(len(mi_list)):
		mid.append([mi_list[i], mi_depth[i]])

	# sort from highest count to lowest count
	mid.sort(key=lambda x:x[1], reverse=True)
	n = len(mid)
	
	num_mi = 0
	
	pid = 0
	cid = 0
	merged = True
	
	# loop until no more merges occure
	while merged:
		# reset merged flag
		merged = False
		# loop through
		for pid in range(n-1):
			if mid[pid] is None:
				continue
			
			# get count of this potential parent	
			pcount = mid[pid][1]
			
			# loop through remaining MIs
			for cid in range(pid+1, n):
				if mid[cid] is None:
					continue
				
				if mid[cid][1]*2 < pcount:
					# compare sequence
					diff = mi_hamming_dist(mid[pid][0], mid[cid][0])
					if diff <= 1:
						# merge to parent. increment parent count
						mid[pid][1] += mid[cid][1]
						pcount = mid[pid][1]
						mid[cid] = None
						merged = True
	
	# loop through the possibly collapsed list of MIs
	for i in range(len(mid)):
		if mid[i] is None:
			continue
			
		# count any MI with more than 1 read supporting them
		if mid[i] > 1:
			num_mi += 1
	
	return num_mi

def mi_hamming_dist(a, b):
	# pair up the bases in a and b
	u = zip(a, b)
	rres = 0

	for l in u:
		if l[0] != l[1]:
			rres += 1
	
	return rres

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
#	fsplit = tmp[8].split("\"")
#	n = len(fsplit)-1
#	i = 0
#	while i < n:
#		key = re.sub(';','',fsplit[i])
#		grow['attrs'][key.strip()] = fsplit[i+1].strip()
#		i += 2

	# split on the semi-colon. it turns out that you can't always count
	# on all values having quotes but the semi-colons are always there.
	# this should parse fine even there aren't quotes on the values and 
	# if there are spaces in the values
	ll = tmp[8].split(";")
	
	# parse out keys and values C style. we'll iterate over the split strings
	# and pull out the pieces as we go
	for sz in ll:
		i0 = 0
		i = 0
		lsz = list(sz)

		if len(lsz) > 0:

			# skip white space at the front
			while lsz[i0]==" ":
				i0 += 1

			# skip to next white space
			i = i0+1
			while lsz[i] != " ":
				i += 1

			name = "".join(lsz[i0:i])

			# parse out the value
			i0 = i+1
			while lsz[i0]==" ":
				i0 += 1

			if lsz[i0]=="\"":
				i0 += 1

			i = i0+1
			while i < len(lsz):
				if lsz[i]=="\"":
					break
				i += 1

			value = "".join(lsz[i0:i])
			
			grow['attrs'][name] = value

#	for k in ll:
#		k = re.sub("^[\s]+", "", k)
#		k = re.sub("[\s]+$", "", k)
#		if len(k) > 0:
#			rres = re.search("^([^\s]+)", k)
#			if rres:
#				name = rres.group(1)
#				k = re.sub(name, "", k)
#				k = re.sub("^[\s]+", "", k)
#				k = re.sub("^\"", "", k)
#				k = re.sub("\"$", "", k)
#				grow['attrs'][name] = k

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
		if szl[0]=="#":
			continue
		grow = gtf_parseline(szl)
		if grow['type'] == "exon":
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
		gid = annot[tid]['gene_id']
		for e in annot[tid]['exons']:
			# tags will be 'gene_id'
			r = region_init(annot[tid]['rname'], e[0], e[1], annot[tid]['strand'], tag=gid)
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
	tid_no_strand = False

	# scan through hash
	for h in lhash:
		if h in lktable:
			# scan through regions in this bucket
			for r0 in lktable[h]:
				rres = compare_regions(r, r0)
				if rres > 0:
					tid_no_strand = r0[REGION_STRAND] != "+" and r0[REGION_STRAND] != "-"
					sense = r[REGION_STRAND]==r0[REGION_STRAND]
					
					if (target_sense is None) or (target_sense == sense) or tid_no_strand:

						# we have a hit!
						has_hits = True
						ovl_len = region_overlap_length(r, r0)

						if ovl_len[0] == 0:
							error_message("found overlap with length 0!")
							print region_str(r), region_str(r0)
							sys.exit(1)

						if r0[REGION_TAG] not in d_hits:
							try:
								d_hits[r0[REGION_TAG]] = { 'length': ovl_len[0] }
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

def sam_mi_from_qname(aln):
	tmp = aln[SAM_QNAME].split(":")
	return tmp[-1]

def sam_qname(aln):
	tmp = aln[SAM_QNAME].split(":")
	n = len(tmp)
	return ":".join(tmp[0:(n-1)])

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

#parser.add_argument('-p', action="store", default=1, type=int, 
#	help="Number of processes to use for finding hits [1]")

parser.add_argument('-o', type=str, default=None, help="Stub for output file. Default is to use the input BAM file name.")

parser.add_argument("-l", "--min-overlap-ratio", default=0.9, type=float, action="store", 
	help="Minimum overlap as a ratio of the read, or read segment, length [0.9]")

parser.add_argument('-q', '--min-mapq', type=int, default=1, 
	help="Minimum MAPQ for counting alignment [1]")

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
		#sys.exit(main_mp(args))
		sys.exit(main(args))
#	except KeyboardInterrupt:
#		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)