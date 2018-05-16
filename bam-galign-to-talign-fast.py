#!/usr/bin/python
#==============================================================================
# bam-galign-to-talign-fast.py
#
# Shawn Driscoll
# 20170608
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This code parses genome aligned reads and intersects them with a GTF
# annotation. The result is SAM output with one line for each read/transcript
# intersection. The SAM data is going to be very stripped down and will contain
# the read name, the FLAG will hold the read's strand orientation relative 
# to the transcript (not the original genome orientation), the target transcript, 
# the original MAPQ, CIGAR and other fields. The actual read and qualities will 
# be removed and replaced with *. 
#==============================================================================

import sys, argparse, math, re, os
from os.path import isfile, expanduser, isdir
from collections import defaultdict
from time import localtime, time
# used for printing exceptions
#import linecache
import traceback
import subprocess as sp
import hashlib
#from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
#import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")
HBIN = 16000

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

# refFlat fields
RF_GENE = 0
RF_TNAME = 1
RF_CHROM = 2
RF_STRAND = 3
RF_TXSTART = 4
RF_TXEND = 5
RF_CDS_START = 6
RF_CDS_END = 7
RF_EXON_COUNT = 8
RF_EXON_STARTS = 9
RF_EXON_ENDS = 10
		
# FLAGs
FLAG_ALIGNED = 0x1
FLAG_PASSED = 0x2
FLAG_ASSIGNED = 0x4
FLAG_MULTIMAPPED = 0x8
FLAG_MULTI_ASSIGNED = 0x10

#==============================================================================
# main
#==============================================================================

def main(args):

	num_parsed = 0
	num_frags = 0
	num_aligned = 0
	num_assigned = 0
	num_assigned_multi = 0
	num_passq = 0
	genome_orientation = 0
	transcript_orientation = 0
	r = None
	read_len = 0
	length_filter_ratio = args.filter_overlap < 1
	
	readidx = {}
	
	if not args.quiet:	
		message("Parsing refFlat annotation", show_time=False)
		
	t0 = time()
	dannot = parse_refflat(args.rflat)
	if not args.quiet:	
		sys.stderr.write("{} sec\n".format(time() - t0))
	
	if not args.quiet:	
		message("Building lookup table", show_time=False)
	t0 = time()
	lktable = gtf_lookup_table(dannot)
	if not args.quiet:	
		sys.stderr.write("{} sec\n".format(time() - t0))
		
	#
	# parse the gtf and build lookup table
	#
	#message("Loading {}".format(args.rflat))
	#dannot, dgid2tid, dgname2tid = gtf_load(args.rflat)
	#message("building index")
	#lktable = gtf_lookup_table(dannot)
	#message("done")
	
	# setup output file
	if args.o is not None:
		if re.search("\.sam", args.o):
			# ok good
			pass
		else:
			args.o += ".sam"
			
		if not args.quiet:	
			if isfile(args.o):
				message("Existing output file will be overwritten")
		
		fout = open(args.o, "w")
	
	else:
		if not args.quiet:	
			message("Results will be written to stdout")
		fout = sys.stdout

	if not args.quiet:	
		message("processing alignments")
	
	fin = open_alignments(args.bam)
	for szl in fin:
		# skip header lines
		if szl[0]=="@":
			continue
		
		# parse
		num_parsed += 1
		if not args.quiet:	
			if (num_parsed % 1000000) == 0:
				progress_message("parsed {}".format(num_parsed))
			
		aln = sam_parse(szl)
		
		if aln[SAM_QNAME] not in readidx:
			readidx[aln[SAM_QNAME]] = 0
		else:
			# multimap flag
			readidx[aln[SAM_QNAME]] |= FLAG_MULTIMAPPED
		
		if sam_unaligned(aln):
			continue
		
		# aligned flag
		readidx[aln[SAM_QNAME]] |= FLAG_ALIGNED
		
		num_aligned += 1
		
		if aln[SAM_MAPQ] < args.filter_mapq:
			continue
		
		# pass mapq flag
		readidx[aln[SAM_QNAME]] |= FLAG_PASSED
		
		num_passq += 1
		
		# get read orientation to genome
		if sam_reversed(aln):
			genome_orientation = -1
		else:
			genome_orientation = 1
		
		read_len = aln[SAM_SEQ]
		
		#sam_aln_extend(aln)
		# get regions for this alignment
		aln_regions = sam_parse_regions(aln)
		# add read length into the output
		aln.append("RL:i:{}".format(read_len))
		hits = {}
		has_hits = False
		
		for r in aln_regions:
			# check for intersections with annotation
			rhk = region_hash(r)
			for k in rhk:
				if k in lktable:
					for r0 in lktable[k]:
						rres = compare_regions(r, r0)
						if rres > 0:
							has_hits = True
							# we got a hit
							tag, eid = r0[REGION_TAG].split(":")
							if tag not in hits:
								# keep track of total overlap length
								hits[tag] = [0, 0, []]
								
								# set alignment's orientation to the transcript
								if r0[REGION_STRAND] == "+":
									if genome_orientation == 1:
										# same
										hits[tag][1] = 1
									else:
										# opposite
										hits[tag][1] = -1
								elif r0[REGION_STRAND] == "-":
									if genome_orientation == 1:
										# opposite
										hits[tag][1] = -1
									else:
										# same
										hits[tag][1] = 1
							
							ovl = region_overlap_length(r, r0)
							hits[tag][0] += ovl[0]
							hits[tag][2].append(eid)
		
		if has_hits:
			# build new output lines
			
			for tid in hits.keys():
				
				# check overlap
				if length_filter_ratio:
					rat = hits[tid][0]*1.0/read_len
					if rat < args.filter_overlap:
						continue

				else:
					if hits[tid][0] < args.filter_overlap:
						continue
				
				# check strand
				if args.f_stranded and hits[tid][1] < 0:
					continue
				
				if args.r_stranded and hits[tid][1] > 0:
					continue
				
				# assigned flag
				if (readidx[aln[SAM_QNAME]] & FLAG_ASSIGNED):
					readidx[aln[SAM_QNAME]] |= FLAG_MULTI_ASSIGNED
					
				readidx[aln[SAM_QNAME]] |= FLAG_ASSIGNED
				
				# make overlap tag
				otag = "OV:i:{}".format(hits[tid][0])
				etag = "EX:Z:{}".format("x".join(map(str, list(set(hits[tid][2])))))
				# make flag
				flag = 0
				if hits[tid][1] < 0:
					flag = 16
				
				ahat = list(aln)
				ahat[SAM_FLAG] = flag
				ahat[SAM_RNAME] = tid
				ahat[SAM_SEQ] = "*"
				ahat[SAM_QUAL] = "*"
				ahat.append(otag)
				ahat.append(etag)
				
				sz = "\t".join(map(str, ahat))
				fout.write(sz + "\n")
				
		
		
	fin.close()
	fout.close()

	if not args.quiet:	
		progress_message("parsed {}".format(num_parsed), last=True)	
	
	if not args.quiet:	
		message("Finished writing alignments. Calculating alignment stats.")
	
	# COMPUTE STATS FOR READS
	num_reads = 0
	num_aligned = 0
	num_passed = 0
	num_assigned = 0
	num_multimapped = 0
	num_multiassigned = 0
	
	for rid in readidx.keys():
		num_reads += 1
		flag = readidx[rid]
		
		if flag & FLAG_ALIGNED:
			num_aligned += 1
		
		if flag & FLAG_PASSED:
			num_passed += 1
		
		if flag & FLAG_MULTIMAPPED:
			num_multimapped += 1
		
		if flag & FLAG_ASSIGNED:
			num_assigned += 1
			
		if flag & FLAG_MULTI_ASSIGNED:
			num_multiassigned += 1
			
	if not args.quiet:	
	
		sys.stderr.write("Run Stats:\n")
		sys.stderr.write("Total Reads:                          {}\n".format(num_reads))
		sys.stderr.write("Total Aligned:                        {} - {:0.2f}% of total\n".format(num_aligned, num_aligned*100.0/num_reads))
		sys.stderr.write("Total Passed MAPQ:                    {} - {:0.2f}% of aligned\n".format(num_passed, num_passed*100.0/num_aligned))
		sys.stderr.write("Total Multi-mapped:                   {} - {:0.2f}% of aligned\n".format(num_multimapped, num_multimapped*100.0/num_aligned))
		sys.stderr.write("Total Assigned:                       {} - {:0.2f}% of passed\n".format(num_assigned, num_assigned*100.0/num_passed))
		sys.stderr.write("Total Assigned to Multiple Features:  {} - {:0.2f}% of assigned\n".format(num_multiassigned, num_multiassigned*100.0/num_assigned))
		
		message("Have a nice day.")
	
	return 0

#==============================================================================
# general functions
#==============================================================================


def open_alignments(fname):
	cmd = "samtools view"
	if re.search("\.sam$", fname):
		cmd += " -S"
		
	cmd += " {}".format(fname)
	
	p = sp.Popen(cmd.split(), stdout=sp.PIPE)
	return p.stdout

#
# refFlat related functions
#

def parse_refflat(fname):
	
	dannot = {}
	
	with open(fname, "r") as fin:
		for szl in fin:
			if szl[0] == "#":
				continue
			
			aln = szl.strip().split("\t")
			tid = aln[RF_TNAME]
			
			dannot[tid] = {
				'gene_name': aln[0],
				'rname': aln[RF_CHROM], 
				'strand': aln[RF_STRAND],
				'start': int(aln[RF_TXSTART]), 
				'end': int(aln[RF_TXEND]), 
				'cds_start': int(aln[RF_CDS_START]),
				'cds_end': int(aln[RF_CDS_END]), 
				'num_exons': int(aln[RF_EXON_COUNT]), 
				'length': 0, 
				'exons': [] 
			}
			
			
			estarts = map(int, re.sub(",$", "", aln[RF_EXON_STARTS]).split(","))
			eends = map(int, re.sub(",$", "", aln[RF_EXON_ENDS]).split(","))
			
			if len(estarts) != len(eends):
				error_message("feature {} has different counts of exon starts and ends".format(aln[RF_TNAME]))
			else:
				for i in range(len(estarts)):
					dannot[tid]['exons'].append([estarts[i], eends[i]])
					dannot[tid]['length'] += eends[i] - estarts[i] + 1
		
		#print tid, dannot[tid]['exons']
	
	return dannot
			

#==============================================================================
# GTF related functions
#==============================================================================

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


parser = argparse.ArgumentParser(description="Rewrites genome alignments as transcriptome alignments. Makes no effort to figure out accurate position of alignments.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser.add_argument('gtf', nargs="+", metavar="wig", type=str, 
#	help="Input wig or wig files to process")

parser.add_argument('rflat', type=str, 
	help="refFlat format annotation.")
parser.add_argument('bam', type=str, 
	help="BAM or SAM input. Must be readable by samtools.")

parser.add_argument('-o', type=str, default=None, 
	help="Output file name. If '.sam' is omitted then it will be added.")

parser.add_argument('--filter-mapq', type=int, default=0, 
	help="Minimum MAPQ for genome alignments.")

parser.add_argument('--filter-overlap', type=float, default=1, 
	help="Minimum overlap between alignment and transcript. If >= 1 then it is interpreted as a number of bases. If 0 < x < 1 then it is interpreted as a ratio of the read length")

parser.add_argument('--quiet', action="store_const", const=True, default=False, 
	help="Do not echo anything during runtime")

strand_group = parser.add_mutually_exclusive_group(required=False)
strand_group.add_argument('--f-stranded', action="store_const", const=True, default=False, 
	help="Keep only alignments that have the same strand orientation as the target transcript")
strand_group.add_argument('--r-stranded', action="store_const", const=True, default=False, 
	help="Keep only alignments that have the opposite strand orientation as the target transcript")

args = parser.parse_args()

# 
# confirm input files exist
#
if not isfile(args.bam):
	error_message("input file {} does not exist".format(args.bam))
	sys.exit(1)

if not isfile(args.rflat):
	error_message("Annotation file {} does not exist".format(args.rflat))
	sys.exit(1)

if __name__ == "__main__":

	try:
		sys.exit(main(args))
#	except KeyboardInterrupt:
#		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)
