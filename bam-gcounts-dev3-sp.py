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
# to features in a GTF. 
#==============================================================================

import sys, argparse, re, os
from os.path import isfile, expanduser
from collections import defaultdict
import hashlib
import subprocess as sp
import numpy as np
from multiprocessing import cpu_count

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")
HBIN = 16000

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	sam_input = False
	output_file = ""

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

	##
	# load the GTF
	message("loading {}".format(args.gtf))
	d_annot, d_gid2tid, d_gn2tid = gtf_load(args.gtf)
	message("found {} transcripts".format(len(d_annot.keys())))
	message("building lookup table")
	d_lktable = gtf_lookup_table(d_annot)
	message("done.")

	for f in args.bam_list:
		alignments = f
				
		sys.stderr.write("\n")
		message("Processing {}".format(alignments))
		
		if not isfile(alignments):
			error_message("cannot find alignment file ({})".format(alignments))
			continue

		if re.search("\.sam$", alignments):
			# sam input
			message("Detected SAM input instead of BAM")
			sam_input = True
	
		## 
		# create output file name
		output_file = re.sub("\.[bams]+$", ".ghits", alignments)
		if isfile(output_file):
			warning_message("Output file exists. It will be overwritten. ({})".format(output_file))
	
		#
		# now we are ready to go!
		results = core(args, alignments, output_file, sam_input, d_annot, d_lktable)


	return 0

#
# runs intersection of alignments with annotation
def core(args, alignments, output_file, sam, annot, lktable):

	hit_table = defaultdict(float)
	frag_mean0 = frag_mean = 0
	frag_n = 0
	max_threads = min([4, cpu_count()/2])

	num_parsed = num_frags = num_aligned = num_assigned = num_unique = num_multi = num_passq = 0
	
	for tid in annot.keys():
		hit_table[tid] = 0.0

	#
	# get the alignments file open and move it
	if alignments == "-":
		# input is sam on stdin
		fin = sys.stdin
	else:

		if sam:
			# input is sam. open it as a file
			fin = open(alignments, "r")

		else:
			# input is a bam. do we need to sort it or can we just load it up?
			
			fifo_stub = hashlib.md5(alignments).hexdigest()

			if args.name_sort_bam:
				# name sort the bam.  we'll fire up a sub process that calls
				# samtools sort and sort into a FIFO then open the FIFO as the
				# input to the main loop
				message("Name sorting the alignments...")
				# make fifo and sort to it in a background process
				if isfile("{}.bam".format(fifo_stub)):
					os.unlink("{}.bam".format(fifo_stub))

				sys.stderr.write("+mkfifo {}.bam\n".format(fifo_stub))
				os.mkfifo("{}.bam".format(fifo_stub))

				# launch sorting process
				#cmd = "samtools sort -n -@ {} -o {} {}".format(max_threads, "{}.bam".format(fifo_stub), alignments)
				cmd = "samtools sort -@ {} -n {} {}".format(max_threads, alignments, fifo_stub)
				sys.stderr.write("CMD: {}\n".format(cmd))
				p0 = sp.Popen(cmd.split())

				# open fifo for reading...
				p1 = sp.Popen("samtools view {}.bam".format(fifo_stub).split(), stdout=sp.PIPE)
				fin = p1.stdout

			else:
				# bam doesn't need to be sorted

				# input is bam, use samtools to read it
				p1 = sp.Popen("samtools view {}".format(alignments).split(), stdout=sp.PIPE)
				fin = p1.stdout

	#
	# alignments are coming in on 'fin'
	#
	
	message("reading alignments")
	
	rname = rname_last = ""
	read_buffer = []


	read_aligned = False
	read_pass_mapq = False

	for szl in fin:

		if szl[0] == "@":
			# header line, skip it
			continue

		aln = samaln_init(szl)
		rname = aln['qname']


		num_parsed += 1

		if (num_parsed % 1e6) == 0:
			message("parsed {} fragments".format(num_frags))
			
		if rname != rname_last and rname_last != "":
			num_frags += 1

			if read_aligned:
				num_aligned += 1
			if read_pass_mapq:
				num_passq += 1

			# is anything aligned in here?
			if len(read_buffer) > 0:
				# we have aligned reads
				rres, frag_len = process_read(read_buffer, annot, hit_table)
				if rres != 0:
					num_assigned += 1

				if rres == 1:
					num_unique += 1
				elif rres == 2:
					num_multi += 1
				
				if frag_len > 0:
					frag_n += 1
					frag_mean0 = frag_mean
					frag_mean = update_mean(frag_mean0, frag_len, frag_n)

			read_buffer = []
			read_aligned = False
			read_pass_mapq = False

		if not samaln_unaligned(aln):
			read_aligned = True		
			if aln['mapq'] >= args.min_mapq:
				read_pass_mapq = True
				# look up hits
				process_alignment(args, lktable, aln)
				read_buffer.append(aln)

		rname_last = rname

	fin.close()

	if len(read_buffer) > 0:
		num_frags += 1

		if read_aligned:
			num_aligned += 1
		if read_pass_mapq:
			num_passq += 1

		# is anything aligned in here?
		if len(read_buffer) > 0:
			num_aligned += 1
			# we have aligned reads
			rres, frag_len = process_read(read_buffer, annot, hit_table)
			if rres != 0:
				num_assigned += 1

			if rres == 1:
				num_unique += 1
			elif rres == 2:
				num_multi += 1

#	message("parsed {} fragments".format(num_frags))
	message("detected fragment mean: {:0.1f}".format(frag_mean))
	
#	print num_parsed, num_frags, num_aligned, num_assigned, num_unique, num_multi

	# print summary
	sys.stderr.write("\n")
	sys.stderr.write("Parsed {} fragments of which:\n".format(num_frags))
	sys.stderr.write("    {} ({:0.1f}%) were aligned and\n".format(num_aligned, num_aligned*100.0/num_frags))
	sys.stderr.write("    {} ({:0.1f}%) had MAPQ >= {}\n".format(num_passq, num_passq*100.0/num_frags, args.min_mapq))
	sys.stderr.write("\n")
	sys.stderr.write("Of {} fragments that passed filters:\n".format(num_passq))
	sys.stderr.write("    {} ({:0.1f}%) were assigned to features with:\n".format(num_assigned, num_assigned*100.0/num_passq))
	sys.stderr.write("      {} ({:0.1f}%) uniquely assigned to genes and\n".format(num_unique, num_unique*100.0/num_assigned))
	sys.stderr.write("      {} ({:0.1f}%) ambiguously assigned to genes\n".format(num_multi, num_multi*100.0/num_assigned))
	sys.stderr.write("\n")	

	message("writing per-transcript hits to {}".format(output_file))
	fout = open(output_file, "w")
	outlines = lout = []
	total_hits = total_fpkm = 0
	
	for tid in sorted(hit_table.keys()):
		
		eff_length = annot[tid]['length']-frag_mean
		hits = hit_table[tid]
		
		if eff_length < 1:
			eff_length = 0
			hits = 0
			eff_hits = 0
		else:
			eff_hits = hits*annot[tid]['length']*1.0/eff_length
		
		total_hits += eff_hits
			
		lout = [
			tid,
		 	annot[tid]['gene_id'], 
		 	annot[tid]['gene_name'],
		 	annot[tid]['rname'], 
		 	annot[tid]['strand'], 
		 	annot[tid]['length'],
		 	eff_length,
		 	hits, 
		 	eff_hits
		]
		outlines.append(lout)
	
	
	for lout in outlines:
		
		if lout[-1] > 0:
			# calculate fpkm
			f = lout[-1]*1e9/(lout[6]*total_hits)
			total_fpkm += f
			lout.append(f)
		else:
			lout.append(0)

	# write header
	fout.write("transcript_id\tgene_id\tgene_name\tchrom\tstrand\tlength\teff_length\thits\teff_hits\tfpkm\ttpm\n")
	for lout in outlines:
		lout.append(lout[-1]*1e6/total_fpkm)
		fout.write("\t".join(map(str, lout)) + "\n")

	fout.close()
	
	return 



#==============================================================================
#
# general functions
#
#==============================================================================

def update_mean(mu0, x, n):
	mu = mu0 + (x-mu0)*1.0/n
	return mu

#
# rbuffer is a set of alignments all with the same name. this function 
# checks hits for the alignments, establishes the alignment weight and 
# assigns the hits to the elements in the passed results dict
def process_read(rbuffer, annot, hit_table):

	pe = False
	result = 0
	hits = set()
	gene_hits = set()

	aln1 = aln2 = []
	frag_len = []

	# sort the reads into left and right
	for aln in rbuffer:
		if samaln_paired(aln):
			pe = True
			if samaln_is_first_mate(aln):
				aln1.append(aln)
			else:
				aln2.append(aln)
		else:
			aln1.append(aln)

	if pe:
		# pair the alignments
		pares = []

		for aln in aln1:
			for mate in aln2:
				if samaln_is_mate(aln, mate):
					# pair the hits
					for hit1 in aln['hits']:
						for hit2 in mate['hits']:
							if hit1[0]==hit2[0]:
								hits.update([hit1[0]])
								# can we use this pair to estimate insert size?
								e1 = list(hit1[1])
								e2 = list(hit2[1])
								if len(e1)==1 and len(e2)==1 and (e1[0]==e2[0]):
									# yes! they both hit the same exon so we can use this one
									# nice and easy
									frag_len.append(abs(aln['tlen']))
							 
					#pares.append([aln, mate])
					

		# paired. now we can make a single set of transcript ids that were hit

#		for pair in pares:
#			t0 = [x[0] for x in pair[0]['hits']]
#			t1 = [x[0] for x in pair[1]['hits']]
#			hits.update(list(set(t0).intersection(set(t1))))

	else:
		frag_len.append(aln['read_len'])
		for aln in aln1:
			t0 = [x[0] for x in aln['hits']]
			hits.update(list(set(t0)))

	# find genes hit
	hits = list(hits)

	if len(hits)==0:
		return result, -1

	for tid in hits:
		gene_hits.update([annot[tid]['gene_id']])

	gene_hits = list(gene_hits)
	num_hits = len(hits)
	num_genes = len(gene_hits)

	result = 1
	if num_genes > 1:
		result = 2

	# weight for assignment of this read to features
	wi = 1.0/num_hits * 1.0/num_genes

	for tid in hits:
		hit_table[tid] += wi

	if len(frag_len) > 0:
		frag_mean = np.mean(frag_len)
	else:
		frag_mean = -1

	return result, frag_mean




def process_alignment(args, lktable, aln):

	# list for hits to this alignment
	aln['hits'] = []
	tmp = []
	dfinal = {}
	aln_strand = "-" if samaln_is_reversed(aln) else "+"
	min_ratio = args.min_overlap_ratio

	target_sense = None
	
	if args.fr_stranded:
		if samaln_paired(aln) and samaln_is_first_mate(aln):
			target_sense = False
		elif samaln_paired(aln) and samaln_is_second_mate(aln):
			target_sense = True
		elif not samaln_paired(aln):
			target_sense = False
	elif args.rf_stranded:
		if samaln_paired(aln) and samaln_is_first_mate(aln):
			target_sense = True
		elif samaln_paired(aln) and samaln_is_second_mate(aln):
			target_sense = False
		elif not samaln_paired(aln):
			target_sense = True
		
	#
	# for each boundary in list, l, make a region and look up hits
	for b in aln['bounds']:
		r = region_init(aln['rname'], b[0], b[1], aln_strand)
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
		rat = dfinal[tid]['length']*1.0/aln['aln_len']
		if rat >= min_ratio:
			# append hit. transcript id and the exons that were hit
			aln['hits'].append([tid, dfinal[tid]['index']])

	return 0


def error_message(sz):
	sys.stderr.write("Error: {}\n".format(sz))

def warning_message(sz):
	sys.stderr.write("Warning: {}\n".format(sz))

def message(sz):
	sys.stderr.write("[bam-gcounts] {}\n".format(sz))

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
					sense = r['strand']==r0['strand']
					
					if (target_sense is None) or (target_sense == sense):

						# we have a hit!
						has_hits = True
						ovl_len = region_overlap_length(r, r0)
				
						if ovl_len[0] == 0:
							error_message("found overlap with length 0!")
							print region_str(r), region_str(r0)
							sys.exit(1)
	
						if r0['tag'] not in d_hits:
							try:
								d_hits[r0['tag']] = { 'index': set([r0['index']]), 'length': ovl_len[0] }
							except:
								print r0['tag'], r0['index'], ovl_len
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
	r = { 'rname':rname, 'start':start, 'end':end, 'strand':strand }
	
	if tag is not None:
		r['tag'] = tag

	if index is not None:
		r['index'] = index

	return r


#
# assuming 'r' is a dict with 'rname', 'start' and 'end' fields
def region_hash(r):
	bin0 = binN = 0

	bin0 = int(r['start'])/HBIN
	binN = int(r['end'])/HBIN

	hout = ["{}:{}".format(r['rname'], bin0)]
	
	if binN > bin0:
		while bin0 < binN:
			bin0 += 1
			hout.append("{}:{}".format(r['rname'], bin0))

	return hout

#
# compare two regions. 
# return value:
# 0 for no overlap
# 1 for overlap
# 2 for identical
def compare_regions(r1, r2):

	# check ref names. if not equal then we're done
	if r1['rname'] != r2['rname']:
		return 0

	# ref names must be equal
	if r1['start']==r2['start'] and r1['end']==r2['end']:
		# starts and ends are identical
		return 2

	# now check for overlap
	if r1['end'] >= r2['start'] and r2['end'] >= r1['start']:
		# overlap!
		return 1

def region_str(r):
	sz_pos = "{}:{}-{}".format(r['rname'], r['start'], r['end'])

	if "tag" in r:
		sz_pos += "|{}".format(r['tag'])

	if "index" in r:
		sz_pos += "|{}".format(r['index'])

	return sz_pos


def region_length(r):
	return r['end']-r['start']+1

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
	len_C = r1['end']-r2['start']+1
	len_D = r2['end']-r1['start']+1

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

#
# parse a SAM record into a dict. if the read is aligned then the boundaries
# of the alignment are calculated as well as the aligned length
def samaln_init(sz):
	aln = {}

	tmp = sz.strip().split("\t")

	aln['qname'] = tmp[0]
	aln['flag'] = int(tmp[1])
	aln['rname'] = tmp[2]
	aln['pos'] = int(tmp[3])
	aln['mapq'] = int(tmp[4])
	aln['cigar'] = tmp[5]
	aln['rnext'] = tmp[6]
	aln['pnext'] = int(tmp[7])
	aln['tlen'] = int(tmp[8])
	aln['seq'] = tmp[9]
	aln['qual'] = tmp[10]

	aln['read_len'] = len(tmp[9])

	# find boundaries of the alignment as well as the aligned length which 
	# would be different from the read length in the event of soft-clipping
	if not samaln_unaligned(aln):
		samaln_calc_bounds(aln)

	aln['attr'] = {}

	# parse out attributes
	if len(tmp) > 11:
		for i in range(11, len(tmp)):
			tparts = tmp[i].split(":")
			value = tparts[2]
			if tparts[1] == "i":
				value = int(value)
			elif tparts[1] == "f":
				value = float(value)

			aln['attr'][tparts[0]] = value

	return aln

#
# return True if aln is the aligned mate of aln0
def samaln_is_mate(aln0, aln):

	if aln0['rnext'] == aln['rnext'] or aln0['rname'] == aln['rname']:
		# usually these are set to '=' when both alignments are in the same 
		# reference
		if aln0['pnext'] == aln['pos'] and aln0['pos'] == aln['pnext']:
			# positions match up
			return True

	return False




#
# functions to check the status of the alignment flag
def samaln_paired(aln):
	return True if (aln['flag'] & 0x1) != 0 else False

def samaln_properly_paired(aln):
	return True if (aln['flag'] & 0x1) != 0 else False

def samaln_unaligned(aln):
	return True if (aln['flag'] & 0x4) != 0 else False

def samaln_mate_unaligned(aln):
	return True if (aln['flag'] & 0x8) != 0 else False

def samaln_is_reversed(aln):
	return True if (aln['flag'] & 0x10) != 0 else False

def samaln_mate_reversed(aln):
	return True if (aln['flag'] & 0x20) != 0 else False

def samaln_is_first_mate(aln):
	return True if (aln['flag'] & 0x40) != 0 else False

def samaln_is_second_mate(aln):
	return True if (aln['flag'] & 0x80) != 0 else False

def samaln_is_secondary(aln):
	return True if (aln['flag'] & 0x100) != 0 else False

#
# compute the boundaries of an alignment
def samaln_calc_bounds(aln):
	aln['bounds'] = []

	left = right = aln['pos']

	op_len = re.split("[A-Z]", aln['cigar'])
	op_typ = re.split("[0-9]+", aln['cigar'])
	
	# drop the extra bit that comes along with the split
	op_len = map(int, op_len[0:(len(op_len)-1)])
	op_typ = op_typ[1:len(op_typ)]
	nop = len(op_typ)

	# loop through
	aln['aln_len'] = 0
	for i in range(nop):
		if op_typ[i] == "M" or op_typ[i] == "D":
			right += op_len[i]
			if op_typ[i] == "M":
				aln['aln_len'] += op_len[i]

		if op_typ[i] == "N":
			# start a new feature and pass the current on
			aln['bounds'].append([left, right-1])
			left = right + op_len[i]
			right = left

	# append feature
	aln['bounds'].append([left, right-1])

	return 0



#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Quantify hits to GTF transcripts from genome alignments.")

parser.add_argument('gtf', type=str, help="GTF annotation corresponding to the reference")

# add these if we want to add multiple bams. maybe we can process them in parallel with 
# processes?

# metavar='bam', nargs='+', 
parser.add_argument('bam_list', type=str, metavar="bam_list", nargs="+", 
	help="Alignments in BAM/SAM format or - to read from stdin. If reading from stdin then input is expected to be SAM")

parser.add_argument("-n", "--name-sort-bam", action="store_const", const=True, default=False,
	help="Name sort the BAM first before running quantification [off]")

#parser.add_argument("-S", "--sam", action="store_const", const=True, default=False, 
#	help="Input is SAM format [off]")

parser.add_argument("-l", "--min-overlap-ratio", default=0.96, type=float, action="store", 
	help="Minimum overlap as a ratio of the read, or read segment, length [0.96]")

parser.add_argument('-q', '--min-mapq', type=int, default=1, 
	help="Minimum MAPQ for counting alignment [1]")

parser.add_argument('--fr-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a reverse stranded library [off]")

parser.add_argument('--rf-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a forward stranded library [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

