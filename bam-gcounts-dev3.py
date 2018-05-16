#!/usr/bin/python
#==============================================================================
# bam-gcounts-dev3.py
#
# Shawn Driscoll
# 20170511
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Clean procedural implemetation for assigning hits from genome alignments
# to features in a GTF. This implementation uses 'multiprocessing' to 
# enable parallel processing of multiple alignment files. 
#==============================================================================

import sys, argparse, re, os
from os.path import isfile, expanduser
from collections import defaultdict
import hashlib
import subprocess as sp
import numpy as np
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock
from math import sqrt
from time import localtime
import random

#==============================================================================
# globals
#==============================================================================

LOG_FILE = "bam_gcounts.run.log"
HOME = expanduser("~")
HBIN = 16000
CHUNK_SIZE_PER_PROC = 100000

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


#==============================================================================
# main
#==============================================================================


def main(args):

	# variables
	sam_input = False
	output_file = ""
	tmpname = ""
	forks = []
	# set number of threads for sorting alignments
	sort_threads = max([cpu_count()/2-args.p, 1])

	# multi-process variables
	tasks = JoinableQueue()
	results = Queue()
	lock = Lock()
	pool = []
	
	output_file_name = ""
	
	# stats

	num_parsed = 0
	num_frags = 0
	num_aligned = 0
	num_passq = 0
	num_assigned = 0
	num_assigned_unique = 0
	num_assigned_multi = 0
	frag_len_mean = 0
	frag_len_stdev = 0
	

	##
	# check stuff
	##
	
	if not args.name_sort_bam:
		sys.stderr.write("\n")
		message("IMPORTANT: If alignments aren't name sorted then multi-mappers won't be handled")
	
	if not isfile(args.gtf):
		error_message("cannot find GTF file ({})".format(args.gtf))
		return 1
		
	if args.fr_stranded and args.rf_stranded:
		error_message("You cannot set both --fr-stranded and --rf-stranded!")
		return 1

	if not isfile(args.alignments):
		error_message("Input alignments file does not exist {}".format(args.alignments))
		return 1

	# setup output file name
	if args.o is not None:
		output_file_name = "{}.ghits".format(args.o)
	else:
		output_file_name = re.sub("\.bam$", ".ghits", args.alignments)
		

	#
	# announce what we're doing
	#
	message("Processing {}".format(args.alignments))
	
	#
	# load the gtf annotation
	#
	message("loading {}".format(args.gtf))
	d_annot, d_gid2tid, d_gn2tid = gtf_load(args.gtf)
	message("found {} transcripts".format(len(d_annot.keys())))

	#
	# deal with alignments. alignments will be buffered by read name and then 
	# written out to sam files for parallel processing
	#
	
	tmpname = hashlib.md5(args.alignments).hexdigest()
	fifo_input = "{}.input.sam".format(tmpname)
	fifo_sort = "{}.sorted.bam".format(tmpname)
	
	fin = None
	
	if args.alignments=="-":
		# input is sam on stdin
		fin = sys.stdin
	else:
		# input could be sam or bam
		if re.search("\.sam$", args.alignments):
			# input is sam. just open it up
			fin = open(args.alignments, "r")
		else:
			# input is bam. does it need to be sorted?
			# either way we need a fifo
			
			# make input fifo
			try:
				os.unlink(fifo_input)
			except:
				pass
			
			if args.name_sort_bam:
				# now we have to fork out to sort. remove existing fifo, if exists
				# and make a new one
				
				try:
					os.unlink(fifo_sort)
				except:
					pass
				
				message("Sorting input alignments by read name")
				pid = os.fork()
				if pid==0:
					os.mkfifo(fifo_sort)
					cmd = "samtools sort -n -@ {} -o {} {}".format(sort_threads, fifo_sort, args.alignments)
					os.system(cmd)
					os._exit(0)
				else:
					forks.append(pid)
				
				# open the sorted alignments bam with samtools as a fifo
				pid = os.fork()
				if pid==0:
					os.mkfifo(fifo_input)
					cmd = "samtools view -o {} {}".format(fifo_input, fifo_sort)
					os.system(cmd)
					os._exit(0)
				else:
					forks.append(pid)
			
			else:
				# no need to sort so just open the bam into a sam
				pid = os.fork()
				if pid==0:
					os.mkfifo(fifo_input)
					cmd = "samtools view -o {} {}".format(fifo_input, args.alignments)
					os.system(cmd)
					os._exit(0)
				else:
					forks.append(pid)
			
			# open 'fifo_input'. it may not exist right away so we loop until it does. 
			while True:
				try:
					fin = open(fifo_input, "r")
					break
				except:
					pass
					

	if fin is None:
		error_message("No input stream was created. Bailing out.")
		
		if len(forks) > 0:
			for pid in forks:
				os.waitpid(pid, 0)
				
		return 1

	if True:
	
		#
		# this processing pipeline bundles reads by read name up to CHUNK_SIZE_PER_PROC 
		# total read bundles. then that buffer is written to a file and the file 
		# name is inserted into the 'tasks' queue for processing. the worker 
		# processes can grab those file names and process the files in parallel.
		#
		
		chunk_size = CHUNK_SIZE_PER_PROC
		chunk_counter = 0
		sam_index = 0
		sam_fname = ""
		sam_handle = None
		
		# create child processes
		for i in range(args.p):
			# create children	
			p = Process(target=kek, args=(tasks, results, args, d_annot, lock,))
			p.daemon = True
			p.start()
			pool.append(p)

		# create first buffer file
		sam_fname = "{}.{}.sam".format(tmpname, sam_index)
		sam_handle = open(sam_fname, "w")
	
		rname = rname_last = ""
		rbuffer = []
		
		# parse alignments and send them out to files for parallel processing
		message("Splitting alignments and processing in parallel")
		if args.r < 1:
			message("Subsampling reads at a rate of {}".format(args.r))
			
		for szl in fin:
			# skip header lines
			if szl[0] == "@":
				continue
						
			# alignment line
			num_parsed += 1
			aln = sam_parse(szl)
			rname = aln[SAM_QNAME]
		
			if rname != rname_last and rname_last != "":
				if args.r == 1 or (args.r < 1 and random.random() < args.r):
					# increment parsed fragment count
					num_frags += 1
					chunk_counter += 1
					
					if len(rbuffer) > 0:
						# increment aligned count
						num_aligned += 1
						
						# buffer has alignments. put them out to a file
						sam_handle.write(buffer_to_string(rbuffer) + "\n")
						
						rbuffer = []
					
					if chunk_counter >= chunk_size:
						
						# close buffer, push into tasks queue and start a 
						# new one
						sam_handle.close()
						tasks.put(sam_fname)
						sam_index += 1
						chunk_counter = 0
						# start a new file
						sam_fname = "{}.{}.sam".format(tmpname, sam_index)
						sam_handle = open(sam_fname, "w")
						
			# check aligned status
			if not sam_unaligned(aln):
				# read is aligned. check mapq
				if aln[SAM_MAPQ] >= args.min_mapq:
					# read passes minimum mapq
					if len(rbuffer)==0:
						# count as passed mapq only once per read name
						num_passq += 1
					rbuffer.append(aln)
			
			rname_last = rname
			if (num_parsed % 1e6) == 0:
				message("parsed {} fragments".format(num_frags))
		
		num_frags += 1
		message("parsed {} fragments [fin]".format(num_frags))
		
		# write final buffer out
		if len(rbuffer) > 0:
			num_aligned += 1
			# buffer has alignments. put them out to a file
			sam_handle.write(buffer_to_string(rbuffer))
			rbuffer = []
			
		# close input stream
		fin.close()
		
		# close all of the sam file
		sam_handle.close()
		tasks.put(sam_fname)
			
		message("Waiting for children to finish counting hits")
	
		# put in the stops
		for p in pool:
			tasks.put(None)
		
		# wait for the children to finish
		tasks.join()
		
		# just to be sure...
		for p in pool:
			p.join()
	
				
	
	# put in a poison pill for the results queue
	results.put(None)
	
	# ok now we have to load all of the individual hits files and build the 
	# final hits table
	
	# parse individual hits files
	message("loading hits from all processes")
	
	# create final dict and populate with all transcript ids in the annotation
	dhits = defaultdict(float)
	for tid in d_annot.keys():
		dhits[tid] = 0
		
	lflen = []
	lfstdev = []
	lflen_n = []
	while True:
		fname = results.get()
		if fname is None:
			break
		
		# parse statistics out of the header
		fin = open(fname, "r")
		hh = fin.readline().strip().split(",")
		dh = {}
		for k in hh:
			tmp = k.split(":")
			dh[tmp[0]] = tmp[1]
			
		num_assigned += int(dh['num_assigned'])
		num_assigned_unique += int(dh['num_assigned_unique'])
		num_assigned_multi += int(dh['num_assigned_multi'])
		lflen.append(float(dh['frag']))
		lfstdev.append(float(dh['stdev'])**2)
		lflen_n.append(float(dh['n']))
		
		# parse the hits
		for szl in fin:
			aln = szl.strip().split("\t")
			dhits[aln[0]] += float(aln[1])
		
		# close file
		fin.close()
		# remove it
		os.unlink(fname)

	# generate final fragment length mean
	frag_len_mean = weighted_mean(lflen, lflen_n)
	frag_len_stdev = sqrt(weighted_mean(lfstdev, lflen_n))

	sys.stderr.write("\n")
	sys.stderr.write("Mean fragment length: {:0.1f} +/- {:0.1f}\n\n".format(frag_len_mean, frag_len_stdev))
	sys.stderr.write("Parsed {} fragments of which:\n".format(num_frags))
	sys.stderr.write("    {} ({:0.1f}%) were aligned and\n".format(num_aligned, num_aligned*100.0/num_frags))
	sys.stderr.write("    {} ({:0.1f}%) had MAPQ >= {}\n".format(num_passq, num_passq*100.0/num_frags, args.min_mapq))
	sys.stderr.write("\n")
	sys.stderr.write("Of {} fragments that passed filters:\n".format(num_passq))
	sys.stderr.write("    {} ({:0.1f}%) were assigned to features with:\n".format(num_assigned, num_assigned*100.0/num_passq))
	sys.stderr.write("      {} ({:0.1f}%) uniquely assigned to genes and\n".format(num_assigned_unique, num_assigned_unique*100.0/num_assigned))
	sys.stderr.write("      {} ({:0.1f}%) ambiguously assigned to genes\n".format(num_assigned_multi, num_assigned_multi*100.0/num_assigned))
	sys.stderr.write("\n")	
	sys.stderr.write("\n")	

	# wait for fork
	#os.waitpid(childs[0], 0)

	try:
		os.unlink(fifo_sort)
	except:
		pass
	
	try:
		os.unlink(fifo_input)
	except:
		pass

	#
	# produce output
	#	
	
	dheader = { 'length':5, 'eff_length':6, 'hits':7, 'eff_hits':8, 'fpkm':9 }
	
	total_hits = 0
	total_fpkm = 0
	outlines = []
	# this first loop starts to build each output line and counts up total hits
	for tid in sorted(dhits.keys()):
		#print "{}\t{}".format(tid, dhits[tid])
		lout = [tid]
		hits = dhits[tid]
		eff_hits = hits

		length = d_annot[tid]['length']
		eff_length = d_annot[tid]['length']-frag_len_mean
		if eff_length <= 4:
			eff_length = 0
			hits = 0
			eff_hits = 0
		else:
			eff_hits = float(hits)*d_annot[tid]['length']/eff_length
		
		lout += [d_annot[tid]['gene_id'], d_annot[tid]['gene_name'], d_annot[tid]['rname'], d_annot[tid]['strand']]
		
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
	
	# make tpm and print!
	
	message("Writing results to {}".format(output_file_name))
	with open(output_file_name, "w") as fout:
		fout.write("transcript_id\tgene_id\tgene_name\tchrom\tstrand\tlength\teff_length\thits\teff_hits\tfpkm\ttpm\n")
		
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
		

	# done!
	

	return 0

#
# kek
# process worker. each worker gets bundled reads, by read name, in chunks
# and then intersects those aligned reads with the annotation. when the 
# worker is finished (task == None) then the hits results are saved to a 
# file and that file name is passed into the result_queue. the main 
# function will handle merging the individual worker results into a 
# final table and output.
def kek(task_queue, result_queue, args, annot, lock):

	# get process name
	name = current_process().name
	# output name
	poutname = "{}.hits".format(name)
	
	# build lookup table for this process
	if args.v:
		message_mp("building lookup table", name, lock)
	lktable = gtf_lookup_table(annot)
	
	hit_table = defaultdict(float)
	
	num_processed = 0
	num_assigned = 0
	num_unique = 0
	num_multi = 0
	
	# initial fragment length values
	frag_len_mean = 200
	frag_len_mean0 = 0
	frag_len_stdev = ((frag_len_mean)*0.8)**2
	frag_len_stdev0 = 0
	frag_len_n = 2
	
	# go!
	while True:
		task = task_queue.get()
		if task is None:
			task_queue.task_done()
			break
		
		# the task is the name of a sam file we need to parse 
		# and find hits for alignments
		if args.v:
			message_mp("Processing {}".format(task), name, lock)
		
		fin = open(task, "r")
		for szl in fin:
			# expand buffered reads out to a list
			rbuffer = buffer_from_string(szl)
			
			# check each alignment for hits
			for aln in rbuffer:
				sam_aln_extend(aln)
				aln0 = list(aln)
				process_alignment(args, lktable, aln)
						
			# process the read
			rres, frag_len = process_read(rbuffer, annot, hit_table, frag_len_mean, sqrt(frag_len_stdev/(frag_len_n-1)))
			
			# handle results
			if rres != 0:
				num_assigned += 1

			if rres == 1:
				num_unique += 1
			elif rres == 2:
				num_multi += 1
			
			if frag_len > 0:
				frag_len_n += 1
				frag_len_mean0 = frag_len_mean
				frag_len_stdev0 = frag_len_stdev
				frag_len_mean, frag_len_stdev = update_mean_and_sd(frag_len_mean0, frag_len_stdev0, frag_len, frag_len_n)
		
			num_processed += 1
#			if (num_processed % 1000000) == 0:
#				message_mp("processed {} fragments".format(num_processed), name, lock)
				
		fin.close()
		os.unlink(task)
		task_queue.task_done()
		if args.v:
			message_mp("finished {}".format(task), name, lock)

	# report final number parsed in this process
	frag_len_stdev0 = sqrt(frag_len_stdev*1.0/(frag_len_n-1))
	message_mp("processed {} fragments. frag len: {:0.1f} +/- {:0.1f}".format(num_processed, frag_len_mean, frag_len_stdev0), name, lock)
		
	# write results from this process to a file
	fout = open(poutname, "w")
	# write a header that contains stats collected
	fout.write("num_assigned:{},num_assigned_unique:{},num_assigned_multi:{},frag:{},stdev:{},n:{}\n".format(num_assigned, num_unique, num_multi, frag_len_mean, frag_len_stdev0, frag_len_n))
	for tid in hit_table.keys():
		fout.write("{}\t{}\n".format(tid, hit_table[tid]))
	fout.close()
	result_queue.put(poutname)
	
	return

#==============================================================================
#
# general functions
#
#==============================================================================

def update_mean_and_sd(mu0, sd0, x, n):
	mu = mu0 + (x-mu0)*1.0/n
	sd = sd0 + (x-mu0)*(x-mu)
	return mu, sd

#
# rbuffer is a set of alignments all with the same name. this function 
# checks hits for the alignments, establishes the alignment weight and 
# assigns the hits to the elements in the passed results dict
def process_read(rbuffer, annot, hit_table, frag_mean, frag_stdev):

	pe = False
	result = 0
	hits = set()
	gene_hits = set()

	aln1 = [] 
	aln2 = []
	frag_len = []
	use_frag_len = True

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
	
		
		for i in range(len(aln1)):
			aln = aln1[i]
			aln_se = [sam_soft_start(aln), sam_end_pos(aln)]
			
			for j in range(len(aln2)):
				mate = aln2[j]
				
				if mate is None:
					continue

				mate_se = [sam_soft_start(mate), sam_end_pos(mate)]
					
				if sam_are_mates(aln, mate):
					# pair the hits
					for hit1 in sam_hits(aln):
						for hit2 in sam_hits(mate):								
							if hit1[0]==hit2[0]:

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
								tid = hit1[0]
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
									if e[0] > left_bounds[REGION_END]:
										# if we're still looking and the exon coordinates are 
										# past the alignment coordinates then we're done
										break
									if left_bounds[REGION_START] <= e[1] and left_bounds[REGION_END] >= e[0]:
										found_start = True
										break
									estart += 1
								
								if found_start:
									eend = estart
								
								while eend < n:
									e = annot[tid]['exons'][eend]
									if e[0] > right_bounds[REGION_END]:
										break
									if right_bounds[REGION_START] <= e[1] and right_bounds[REGION_END] >= e[0]:
										found_end = True
										break
									eend += 1
								
								if found_start and found_end:
									# found both start and end so we're good
								
									if eend == estart:
										isize = end-start+1
									else:
										# translate it
										tstart = max([0, start - annot[tid]['exons'][estart][0]]) + annot[tid]['texons'][estart][0]
										tend = max([0, end - annot[tid]['exons'][eend][0]]) + annot[tid]['texons'][eend][0]
										isize = tend-tstart+1
									
									# use insert size to filter out inappropriate (emperically) insert sizes.
									# this helps cut down on mis assignment of reads to transcripts									
									isize_z = (isize-frag_mean)*1.0/frag_stdev
									if isize > 0 and isize_z < 5:
										frag_len.append(isize)
										hits.update([hit1[0]])
																	
								else:
									# one reason for this is if you have paired-end reads and one of the reads
									# overlaps a transcript but the other does not. so like an example of an aligned
									# fragment that doesn't belong to this transcript
									
									pass
									#print "-- ERROR"
									#print start,end,estart,eend
									#print left_bounds, right_bounds
									#print annot[tid]['exons']
									#print annot[tid]['texons']
									#print "--"
									
							
							# paired or not
						# done with mate hits
					# done with aln hits
					# done with paired mates. null out the second mate so it gets skipped over
					# faster next time around
					aln2[j] = None	
		

	else:
		#----------------------------------------------------------------------
		# 
		# single-end condition
		#
		#----------------------------------------------------------------------
		for aln in aln1:
			# use aligned length as fragment length
			frag_len.append(aln[SAM_ALNLEN])
			# make list of the hits
			t0 = [x[0] for x in sam_hits(aln)]
			hits.update(list(set(t0)))

	# find genes hit
	hits = list(hits)

	if len(hits)==0:
		return result, -1

	for tid in hits:
		gene_hits.update([annot[tid]['gene_name']])

	gene_hits = list(gene_hits)
	num_hits = len(hits)
	num_genes = len(gene_hits)

	result = 1
	if num_genes > 1:
		result = 2

	# weight for assignment of this read to features
	wi = 1.0/num_hits * 1.0/num_genes**2

	for tid in hits:
		hit_table[tid] += wi

	if len(frag_len) > 0:
		frag_mean = np.mean(frag_len)
	else:
		frag_mean = -1

	return result, frag_mean



#
# check alignment, aln, for valid hits to the GTF annotation. the 
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

	#
	# merge hits down to a single set of transcripts
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
						# collision. combine them into a single hit
						dfinal[tid]['index'].update(d[tid]['index'])
						dfinal[tid]['length'] += d[tid]['length']

	#
	# evaluate hit overlap length relative to acceptable ratio
	for tid in dfinal.keys():
		rat = dfinal[tid]['length']*1.0/aln[SAM_ALNLEN]
		if rat >= min_ratio:
			# append hit. transcript id and the exons that were hit
			aln[hit_index].append([tid, dfinal[tid]['index']])

	return 0


def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))


def message_mp(sz, name, lock):
	lock.acquire()
	sys.stderr.write("[{}] {}\n".format(name, sz))
	lock.release()
	return


def weighted_mean(mu, wi):
	
	sum_prod = 0
	sum_wi = 0
	
	for i in range(len(mu)):
		sum_prod += mu[i]*wi[i]
		sum_wi += wi[i]
	
	return sum_prod/sum_wi


def buffer_to_string(buf):
	sz = ""
	nl = len(buf)
	for i in range(nl):
		l = map(str, buf[i])
		n = len(l)
		if i > 0:
			sz += "\t"
		sz += "{}\t{}".format(n, "\t".join(l))

	return sz

def buffer_from_string(sz):
	tmp = sz.strip().split("\t")
	buf = []

	i = 0
	idx = 0
	while i < len(tmp):
		# get number of fields
		n = int(tmp[i])	
		# slice from next to i+1+n
		buf.append(list(tmp[(i+1):(i+1+n)]))
		# convert fields into ints
		for k in [SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_TLEN, SAM_RLEN, SAM_ALNLEN]:
			buf[idx][k] = int(buf[idx][k])
		idx += 1
		i = i + 1 + n

	return buf

#==============================================================================
#
# GTF related functions
#
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

def gtf_row_tostring(obj):
	lout = [obj['rname'], obj['db'], obj['type'], str(obj['start']), str(obj['end']), obj['strand']]
	return "\t".join(lout)

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

	with open(fname, "r") as fin:
		for szl in fin:
			
			grow = gtf_parseline(szl)
			
			tid = gtf_transcript_id(grow)
			gid = gtf_gene_id(grow)
			gname = gtf_gene_name(grow)
			
			if (tid is not None) and (gid is not None):
					
				if tid not in dannot:

					if gname is None:
						#warning_message("Transcript, {}, is not annotated with a gene name. Using gene id as gene name".format(tid))
						gname = gid

					dannot[tid] = {
						'rname':grow['rname'],
						'strand':grow['strand'],
						'db':grow['db'],
						'gene_id':gid,
						'gene_name':gname,
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
			
			else:
				if tid is None:
					warning_message("Transcript without a transcript id was found: {}".format(gtf_row_tostring(grow)))
				elif gid is None:
					warning_message("Transcript, {}, is not annotated with a gene_id and was skipped".format(tid))

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
			r = region_init(annot[tid]['rname'], e[0], e[1], annot[tid]['strand'], tag=tid, index=eidx)
			eidx += 1
			h = region_hash(r)
			for hid in h:
				# insert the region in each binf
				lktable[hid].append(r)

	return lktable


#
# look up region r in the lookup table. return info from each hit including 
# the tag and index values
def gtf_find_hits(lktable, r, target_sense):

	# hash the region, r
	lhash = region_hash(r)
	# hits list
	d_hits = {}
	has_hits = False

	# scan through hashs
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
#
# SAM related functions
#
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

	for i in [SAM_FLAG, SAM_POS, SAM_PNEXT]:
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

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Quantify hits to GTF transcripts from genome alignments.")

parser.add_argument('gtf', type=str, help="GTF annotation corresponding to the reference")

parser.add_argument('alignments', type=str, 
	help="SAM or BAM alignments. You may also enter '-' to pass alignments via stdin. NOTE: alignments must be name sorted or else your results will be wacked.")

# add these if we want to add multiple bams. maybe we can process them in parallel with 
# processes?

# metavar='bam', nargs='+', 
#parser.add_argument('bam_list', type=str, metavar="bam_list", nargs="+", 
#	help="Alignments in BAM/SAM format or - to read from stdin. If reading from stdin then input is expected to be SAM")

parser.add_argument('-o', action="store", default=None, type=str, 
	help="Stub for output file. If not supplied then the alignment file name is used")

parser.add_argument("-p", action="store", default=1, type=int,
	help="Process this many files in parallel (if you specified more than one file...)")

parser.add_argument("-n", "--name-sort-bam", action="store_const", const=True, default=False,
	help="Name sort the BAM first before running quantification [off]")

parser.add_argument('-r', action="store", default=1, type=float, 
	help="Subsample alignments at this rate. To quantify 10%% of the reads set this to 0.1. [1]")

#parser.add_argument("-S", "--sam", action="store_const", const=True, default=False, 
#	help="Input is SAM format [off]")

parser.add_argument("-l", "--min-overlap-ratio", default=0.96, type=float, action="store", 
	help="Minimum overlap as a ratio of the read, or read segment, length [0.96]")

parser.add_argument('-q', '--min-mapq', type=int, default=1, 
	help="Minimum MAPQ for counting alignment [1]")

#parser.add_argument('--depth-counts', action="store_const", const=True, default=False, 
#	help="Produce depth and average depth instead of read counts.")

#parser.add_argument('--exon-quant', action="store_const", const=True, default=False, 
#	help="Produce exon-level quantification instead of transcript level")

parser.add_argument('--fr-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a reverse stranded library [off]")

parser.add_argument('--rf-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a forward stranded library [off]")

parser.add_argument('-v', action="store_const", const=True, default=False, 
	help="Verbose output [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

