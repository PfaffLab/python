#!/usr/bin/python
#==============================================================================
# bam-genome-expr-quant.py
#
# Shawn Driscoll
# 20180118
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse genome aligned reads and produce TCC output like kallisto pseudo
#==============================================================================

import sys
import argparse
import math
import re
import os
import traceback
from collections import defaultdict
from time import localtime, time
from os.path import isfile
import hashlib
from random import random
from numpy.random import choice
import copy

#==============================================================================
# globals
#==============================================================================

COUNT_FILE = "counts.mex"
EC_FILE = "counts.ec"
SAMPLE_FILE = "counts.samples"

HBIN = 16000

# sam flags
SAM_PAIRED = 0x1
SAM_PROPERLY_PAIRED = 0x2
SAM_UNALIGNED = 0x4
SAM_MATE_UNALIGNED = 0x8
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
SAM_RLEN = 9
SAM_ALNLEN = 10

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

	is_pe = False
	is_bam = False
	tmp_base = ""
	aln_file = ""
	offset = 0
	read_index = None
	ec_idx = 0
	dec = {}
	dec_counts = {}
	file_index = 0
	counts_out = []
	file_stats = {}
	
	file_err = False
	missing_files = []
	
	for fname in args.alignments:
		file_stats[fname] = {'total_reads':0, 'assigned_reads':0}
		if  not isfile(fname):
			file_err = True
			missing_files.append(fname)

	if file_err:
		error_message("One or more input files do not exist!")
		sys.stderr.write(", ".join(missing_files) + "\n")
		return 1	

	if not isfile(args.annotation):
		error_message("Input annotation does not exist")
		return 1

	##
	## load the annotation
	message("Loading annotation")
	## this also returns a pre-populated dict with all annotated transcripts.
	## we can use this as an initial equivalence class table which will be
	## added to as necessary in the quantificaiton loop which follows
	dannot, dtargets = parse_refflat(args.annotation)
	dec = copy.deepcopy(dtargets)
	ec_idx = len(dec.keys())
	
	# build lookup
	message("building lookup table")
	lktable = build_lookup_table(dannot)
	
	##
	## process all alignment files
	for file_index in range(len(args.alignments)):
		
		alignment_file = args.alignments[file_index]
		message("Processing {}".format(alignment_file))
		
		tmp_base = hashlib.md5(alignment_file).hexdigest()
	
		##
		# deal with SAM/BAM situation
		if re.search("\.bam$", alignment_file):
			is_bam = True
			message("Converting BAM alignments to SAM")
			aln_file = "{}.sam".format(tmp_base)
			rres = bam2sam(alignment_file, aln_file)
	
			if not isfile(aln_file):
				error_message("Conversion from BAM to SAM must have failed.")
				return 1
	
		else:
			aln_file = alignment_file
	
		##
		## open alignments and index reads
		##
		message("Indexing reads", show_time=False)
		read_index = build_read_index(aln_file)
		tmp = read_index.keys()
		file_stats[alignment_file]['total_reads'] = len(tmp)
		
		if args.samplerate > 1 or (args.samplerate < 1 and args.samplerate > 0):
			# need to grab a subset of the read_index keys
			
			if args.samplerate < 1:
				N = int(math.floor(len(read_keys) * args.samplerate))
				message("Using random subset of {} reads".format(N), show_time=False)
			else:
				N = int(math.floor(args.samplerate))
				message("Using random subset of {} reads".format(N), show_time=False)
			
			random_idx = choice(len(tmp), N)
			
			read_keys = [tmp[i] for i in random_idx]
		
		else:
			
			read_keys = tmp
				
			
	
		##
		## now we have the read names indexed and we can do stuff
		##
		message("Counting hits to equivalence classes", show_time=False)
		with open(aln_file, "r") as fin:
			
			##
			## initialize the dec_counts dict from the current dec dict
			##
			dec_counts = {}
			multi_hits = set()
			for ecid in dec.keys():
				dec_counts[ecid] = 0
			
			rcount = 0
			
			for rname in read_keys:
				# do stuff
				rcount += 1
				if (rcount % 1000000) == 0:
					progress_message("processed {} reads".format(rcount))
					
				offsets = read_index[rname]
				mate1 = []
				mate2 = []
	
				for offset in offsets:
	
					if args.samplerate < 1:
						rx = random()
						if rx > args.samplerate:
							continue
	
					fin.seek(offset)
					aln = sam_parse_line(fin.readline())
					
					if (aln[SAM_FLAG] & SAM_UNALIGNED) == 0:
						if aln[SAM_FLAG] & SAM_FIRST_MATE:
							mate1.append(aln)
						elif aln[SAM_FLAG] & SAM_SECOND_MATE:
							mate2.append(aln)
						else:
							# single end data...
							mate1.append(aln)
	
				hits = deal_with_read(mate1, mate2, dannot, lktable, args)
	
				if len(hits) > 0:
					
					file_stats[alignment_file]['assigned_reads'] += 1
					
					ec_id = ",".join(sorted(hits))
					
					if len(hits) > 1:
						multi_hits.add(ec_id)
						
					if ec_id not in dec:
						## add new equivalence class to the dec and dec_counts dicts. 
						## also increment the ec index
						dec[ec_id] = ec_idx
						dec_counts[ec_id] = 0
						ec_idx += 1
	
					dec_counts[ec_id] += 1
			
			progress_message("processed {} reads".format(rcount), last=True)
	
		if is_bam and isfile("{}.sam".format(tmp_base)):
			message("Removing temporary files", show_time=False)
			os.unlink("{}.sam".format(tmp_base))
	
		##
		## note that right here we could resolve the multimappers
		##
		
		multi_hits = list(multi_hits)
		if len(multi_hits) > 0:
			message("resolving multi-mappers", show_time=False)
			
			mm_weights = {}
			mm_did = {}
			mm_dlengths = {}
			mm_reads = 0
			for ecid in multi_hits:
				mm_reads += dec_counts[ecid]
				
				ecid_parts = ecid.split(",")
				# set initial weights
				mm_weights[ecid] = sumnorm([1 for i in range(len(ecid_parts))])
				# keep the list of split up target names
				mm_did[ecid] = ecid_parts
				mm_dlengths[ecid] = [dannot[tid]['length'] for tid in ecid_parts]

			message("{} reads in {} equivalence classes".format(mm_reads, len(multi_hits)), show_time=False)

			# get all of the non-multimap hits to this file
			dhits = {}
			u_reads = 0
			for ecid in dec_counts.keys():
				tmp = ecid.split(",")
				if len(tmp)==1:
					u_reads += dec_counts[ecid]
					dhits[ecid] = dec_counts[ecid]
			
			message("{} reads mapped uniquely. unique + multi = {}".format(u_reads, u_reads+mm_reads), show_time=False)
			
			dhits_base = copy.deepcopy(dhits)
			dhits_hat = copy.deepcopy(dhits)
			
			# dhits will update at each iteration. dhits will be used to 
			# calculate weights
			
			# calc weights from dhits
			# update hits in dhits_hat
			# copy dhits_hat to dhits
			# copy dhits_base to dhits_hat
			# when finished, dhits is the final result
			
			max_iter = 1000
			iter = 0
			tol = 1e-9
			sse = 0 
			sse0 = 0
			
			while True:
				iter += 1
				sse = 0
				
				for ecid in multi_hits:
					# initial weight for this iterations
					w0 = mm_weights[ecid]
					
					# update the weights based on length normalized counts
					# from dhits
					w = list(w0)
					for i in range(len(w)):
						w[i] = dhits[mm_did[ecid][i]] * 1.0 / mm_dlengths[ecid][i]
					w = sumnorm(w)
					
					# use the weights to update counts. if the weight vector is all zeros 
					# it means the individual transcripts don't have any unique hits at all dispite
					# having some shared hits. so we'll distribute the shared hits evenly
					if sum(w)==0:
						w = sumnorm([1 for i in range(len(w0))])
					
					# update weight vector
					mm_weights[ecid] = list(w)
					
					# compute change in weights
					sse += sum([(w0[i]-w[i])**2 for i in range(len(w))]) * 1.0 / len(w)
																	
					for i in range(len(w)):
						# add in weighted proportion the multimap count to the target transcript
						dhits_hat[mm_did[ecid][i]] += w[i]*dec_counts[ecid]
										
					
				# copy
				dhits = copy.deepcopy(dhits_hat)
				dhits_hat = copy.deepcopy(dhits_base)
				
				# check error
				sse = sse * 1.0 / len(multi_hits)
				
				progress_message("iteration: {}; mse: {}".format(iter, sse))
				
				if iter >= max_iter:
					message("\ndid not converge", show_time=False)
					break
				
				if abs(sse0-sse) < tol:
					message("\nconverged", show_time=False)
					break
				
				sse0 = sse


			##
			## write this file's counts to the 'counts_out' list
			message("Reporting hits", show_time=False)
			for ecid in dhits.keys():
				if dhits[ecid] > 0:
					# append row with the ec class index, file index and the count
					counts_out.append([dec[ecid], file_index, dhits[ecid]])
			
		else:
	
			##
			## write this file's counts to the 'counts_out' list
			message("Reporting hits", show_time=False)
			for ecid in dtargets.keys():
				if dec_counts[ecid] > 0:
					# append row with the ec class index, file index and the count
					counts_out.append([dec[ecid], file_index, dec_counts[ecid]])
	
	##
	## finished processing files
	message("Finished processing alignment files. Writing output files.")
	
	with open(COUNT_FILE, "w") as fout:
		for lout in counts_out:
			fout.write("\t".join(map(str, lout)))
			fout.write("\n")
	
	##
	## make an ordered list of the ecs
	ec_out = []
	for ecid in dtargets.keys():
		ec_out.append([dtargets[ecid], ecid])
	
	# sort it
	ec_out.sort(key=lambda x: x[0])
	
	with open(EC_FILE, "w") as fout:
		for lout in ec_out:
			fout.write(lout[1])
			fout.write("\n")
		
	##
	## write samples file
	with open(SAMPLE_FILE, "w") as fout:
		for fname in args.alignments:
			lout = [fname, file_stats[fname]['total_reads'], file_stats[fname]['assigned_reads']]
			fout.write("\t".join(map(str, lout)) + "\n")
	
	
	message("Done.")

	return 0


def sumnorm(v):
	ss = sum(v)
	if ss > 0:
		vhat = [a*1.0/ss for a in v]
	else:
		vhat = list(v)
		
	return(vhat)

def format_index(v):
	
	# explode the number out into a list of characters
	vhat = list(str(v))
	n = len(vhat)
	if n > 3:
		pass
	
	return
	

def deal_with_read(mate1, mate2, annot, lktable, args):

	rres0 = []
	rres = []

	pe = False
	if len(mate1) > 0:
		if (mate1[0][SAM_FLAG] & SAM_PAIRED) != 0:
			pe = True

#	if len(mate1) > 0 and len(mate2) > 0:
	if pe:
		
		rres0 = process_pe(mate1, mate2, annot, lktable, args)
		# check length of overlaps
		if args.f < 1:
			for hit in rres0:
				total_ratio = (hit[1] + hit[2]) / 2.0
				if total_ratio >= args.f:
					rres.append(hit[0])

		else:
			rres = list(rres0)

	elif len(mate1) > 0:
		rres0 = process_se(mate1, annot, lktable, args)		
		if args.f < 1:
			for hit in rres0:
				if hit[1] >= args.f:
					rres.append(hit[0])

		else:
			rres = list(rres0)

	# return list of unique targets that have passed filters
	return list(set(rres))


def process_pe(mate1, mate2, annot, lktable, args):

	lpares = []
	total_hits = []
	fstrand = args.f_stranded
	
	for i in range(len(mate1)):
		aln1 = mate1[i]

		if aln1 is None:
			continue
		
		for j in range(len(mate2)):
			aln2 = mate2[j]

			if aln2 is None:
				continue

			if sam_are_mates(aln1, aln2):
				# good!
				lpares.append([aln1, aln2])

				mate1[i] = None
				mate2[j] = None
				continue

	## mates are paired now so we can collect all hits
	left_hits = []
	left_hits0 = []
	right_hits = []
	right_hits0 = []

	for pp in lpares:
		
		# check mapq
		if pp[0][SAM_MAPQ] < args.q or pp[1][SAM_MAPQ] < args.q:
			continue
		
		left_regions = sam_parse_regions(pp[0])
		left_hits0 = []
		for r0 in left_regions:
			left_hits0 += lktable_search(r0, lktable)

		right_regions = sam_parse_regions(pp[1])
		right_hits0 = []
		for r0 in right_regions:
			right_hits0 += lktable_search(r0, lktable)

		left_hits = collapse_hits(left_hits0, pp[0][SAM_ALNLEN])
		right_hits = collapse_hits(right_hits0, pp[1][SAM_ALNLEN])

		##
		## collapse the hits for this pair to only targets that they have in common
		tmp = collapse_pe_hits(left_hits, right_hits)
		
		##
		## deal with strand of this alignment
		if args.f_stranded or args.r_stranded:
			# have to check strand of the alignments relative to the strand of 
			# the transcripts
			
			# get strand of first mate
			astrand = "+"
			if pp[0][SAM_FLAG] & SAM_REVERSED:
				astrand = "-"
	
			total_hits += strand_filter_hits(astrand, fstrand, tmp, annot)
		
		else:
			total_hits += tmp
	
	return total_hits


def process_se(mate1, annot, lktable, args):

	total_hits = []
	hits0 = []
	fstrand = args.f_stranded
	
	for i in range(len(mate1)):
		
		# check mapq
		if mate1[i][SAM_MAPQ] < args.q:
			continue
		
		hits0 = []
		aln_regions = sam_parse_regions(mate1[i])
		for r0 in aln_regions:
			hits0 += lktable_search(r0, lktable)
		
		tmp = collapse_hits(hits0, mate1[i][SAM_ALNLEN])
	
		if args.f_stranded or args.r_stranded:
			# check strand of the hits
			# get strand of this alignment
			astrand = "+"
			if mate1[i][SAM_FLAG] & SAM_REVERSED:
				astrand = "-"
			
			total_hits += strand_filter_hits(astrand, fstrand, tmp, annot)
		
		else:
			
			total_hits += tmp
			
	


def collapse_hits(lhits, aligned_length):
	# list of 2 item lists that include target and overlap length
	dtargets = defaultdict(int)
	rres = []

	if len(lhits) == 0:
		return rres

	for hit in lhits:
		dtargets[hit[0]] += hit[1]

	for tid in dtargets.keys():
		rres.append([tid, dtargets[tid]*1.0/aligned_length])

	return rres

def collapse_pe_hits(lhits, rhits):

	dlhits = {}
	drhits = {}

	rres = []

	for h in lhits:
		dlhits[h[0]] = h[1]

	for h in rhits:
		drhits[h[0]] = h[1]

	for tid in dlhits.keys():
		if tid in drhits:
			rres.append([tid, dlhits[tid], drhits[tid]])

	return rres

def strand_filter_hits(strand0, forward, lhits, annot):

	rres = []

	for h in lhits:
		hstrand = annot[h[0]]['strand']

		if forward:
			if strand0 == hstrand:
				rres.append(h)

		else:
			if strand0 != hstrand:
				rres.append(h)

	return rres

def build_read_index(sam_file):

	offset = 0
	read_index = defaultdict(list)
	line_idx = 0

	with open(sam_file, "r") as fin:

		for szl in fin:
			if szl[0] != "@":
				line_idx += 1
				if (line_idx % 1000000) == 0:
					progress_message("{} lines".format(line_idx))
				# we have an alignment!
				aln = szl.strip().split("\t")
				qname = aln[0]
				if (int(aln[1]) & 0x1) != 0:
					is_pe = True

				read_index[qname].append(offset)

			offset += len(szl)
	
	progress_message("{} lines".format(line_idx), last=True)
	message("{} distinct reads".format(len(read_index.keys())), show_time=False)

	return read_index

###############################################################################
##
## ANNOTATION STUFF
##
###############################################################################

def parse_refflat(fname):
	
	dannot = {}
	dtid_index = {}
	tid_index = 0
	ltid = []
	
	with open(fname, "r") as fin:
		for szl in fin:
			if szl[0] == "#":
				continue
			
			aln = szl.strip().split("\t")
			tid = aln[RF_TNAME]
			
			dtid_index[tid] = tid_index
			tid_index += 1
			
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


	#
	# pass through the annotation and sort the exons by position. also create 
	# transcript versions of the coordinates. this will be used to calculate
	# insert size for PE alignments
#	for tid in dannot.keys():
#		dannot[tid]['exons'].sort(key=lambda x: x[0])
#		# create transcript coordinate version of the exons
#		dannot[tid]['texons'] = []
#		left = 1
#		right = left
#		for e in dannot[tid]['exons']:
#			l = e[1]-e[0]+1
#			right = left + l - 1
#			dannot[tid]['texons'].append([left, right])
#			left = right+1
			
	return dannot, dtid_index

#
# build a quick-lookup table for exon features from a loaded gtf
def build_lookup_table(annot):

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

def lktable_search(r, lktable):
	# hash region
	rh = region_hash(r)
	hits = []

	for h in rh:
		if h in lktable:
			for r0 in lktable[h]:
				rres = region_overlaps(r, r0)
				if rres > 0:
					# got a hit. append name and length of overlap.
					tid = r0[REGION_TAG].split(":")
					hits.append([tid[0], rres])

	return hits

###############################################################################
##
## REGION CODES
##
###############################################################################


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


def region_length(r):
	return r[REGION_END]-r[REGION_START]+1

def region_overlaps(r1, r2):

	if r1[REGION_RNAME] == r2[REGION_RNAME]:
		if r2[REGION_START] <= r1[REGION_END] and r1[REGION_START] <= r2[REGION_END]:
			return region_overlap_length(r1, r2)

	return -1

def region_overlap_length(r1, r2):
	l1 = region_length(r1)
	l2 = region_length(r2)
	l3 = r1[REGION_END] - r2[REGION_START] + 1
	l4 = r2[REGION_END] - r1[REGION_START] + 1
	rres = min([l1, l2, l3, l4])
	return rres

###############################################################################
##
## SAM FILE CODES
##
###############################################################################

# 
# split a sam line
def sam_parse_line(sz):
	aln = sz.strip().split("\t")

	for i in [SAM_FLAG, SAM_POS, SAM_PNEXT, SAM_MAPQ, SAM_TLEN]:
		aln[i] = int(aln[i])
	
	# replace the sesquence with the length of the sequence
	aln[SAM_SEQ] = len(aln[SAM_SEQ])
	# replace the quals with the aligned length
	aln[SAM_QUAL] = sam_aligned_length(aln)
		
	return aln


# expecting that passed object is the list that was split from a single
# sam alignment line
def sam_aligned_length(laln):
	cigar = laln[SAM_CIGAR]
	r = re.findall("([0-9]+)[MD]", cigar)
	return sum(map(int, r))


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


##
# function to call samtools for conversion of alignments
def bam2sam(infile, outfile):
	cmd = "samtools view -o {} {}".format(outfile, infile)
	rres = runcmd(cmd)
	return rres

###############################################################################
##
## GENERAL STUFF
##
###############################################################################

def runcmd(cmd):
	sys.stderr.write("CMD: {}\n".format(cmd))
	return os.system(cmd)

def progress_message(sz, last=False):
	message_length = len(sz)
	field_width = 80
	
	pad_length = field_width - message_length
	if pad_length > 0:
		tmp = [" " for i in range(pad_length)]
		sz = sz + "".join(tmp)
		 
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {} ".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))

def error_message(sz):
	sys.stderr.write("Error: {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

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

parser = argparse.ArgumentParser(description="Quantify one or more genome alignement files vs a refFlat annotation. Output is similar to kallisto 'pseudo' or cellranger's quantification.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('annotation', type=str, 
	help="Transcript annotation in refFlat format")
parser.add_argument('alignments', type=str, nargs="+", 
	help="Alignments in SAM or BAM format for processeing.")

parser.add_argument('-f', type=float, default=0.25, 
	help="Minimum overlap ratio between read and transcript. If >= 1 then it is interpreted as a number of bases. If 0 < x < 1 then a ratio of the aligned read length.")

parser.add_argument('-q', type=int, default=1, 
	help="Minimum MAPQ for accepting alignments.")

parser.add_argument('--samplerate', type=float, default=1, 
	help="Sample rate for reads. If > 1 then it is taken as a specific number of reads to process (randomly selected). If 0 < x < 1 then it is the ratio of aligned reads to quantify.")

strand_group = parser.add_mutually_exclusive_group(required=False)
strand_group.add_argument('--f-stranded', action="store_const", const=True, default=False, 
	help="Library is first-stranded (either SE or PE). First mates would be the same orientation as the transcripts.")
strand_group.add_argument('--r-stranded', action="store_const", const=True, default=False, 
	help="Library is second-stranded (either SE or PE). First mates would be the opposite orientation of the transcripts. (normal stranded case)")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

