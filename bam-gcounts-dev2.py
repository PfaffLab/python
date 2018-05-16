#!/usr/bin/python
#==============================================================================
# bam-gcounts-dev.py
#
# Shawn Driscoll
# 20170505
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Processing genomes alignments to assign hits to features in a 
# gtf file. this is a 'dev' version but it seems to work really well. 
# This pipeline parses the sam alignments and immediatly looks up intersections
# in the lookup table. it deals with bundling reads and assigning hits afterwards.
# the idea is that this may be parallelizable since we can use a
# consumer/producer setup on bundled and intersected reads 
#==============================================================================

import sys, argparse, math, re, hashlib
from os.path import isfile, expanduser
from os import system
import numpy as np
# use this as a way to combine all dicts
from collections import defaultdict
import subprocess as sp

import multiprocessing as mp
import timeit

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")
HBIN = 16000
STRAND_POS = "+"
STRAND_NEG = "-"

ASSIGNED = 1
ASSIGNED_GENE_UNIQUE = 2
ASSIGNED_GENE_MULTI = 4

MAX_MAP_BUFFER = 1000

#==============================================================================
# main
#==============================================================================	

def pworker(task_queue, result_queue):
	while True:
		item = task_queue.get()
		if item is None:
			task_queue.task_done()
			break
		
		rres = process_read(item)
		result_queue.put(rres)
		task_queue.task_done()
	
	return

def main(args):
	
	# variables
	annot = lktable = None
	lset = []
	aln = []
	num_parsed = 0
		
	strand_mode = None
	if args.fr_stranded:
		strand_mode = "FR"
	elif args.rf_stranded:
		strand_mode = "RF"
	
	aln_buffer = []
	map_buffer = []
	
	rname = None
	rname_last = None
	
	fifo_stub = None
	fout_name = None
	
	pid = None
	childs = []

	total_searches = 0
	total_search_time = 0

	hit_factory = None

	# -------------------------------------------------------------------------
	# load GTF annotation
	# -------------------------------------------------------------------------

	sys.stderr.write("loading GTF annotation...\n")
	annot, lktable = parse_gtf(args.gtf)
	sys.stderr.write("found {} features\n".format(len(annot.keys())))
	
	# -------------------------------------------------------------------------
	# name sort the bam or not...
	# -------------------------------------------------------------------------
	
	# loop through BAM files specified at the command line. this saves time 
	# for multiple files since we only have to load the GTF once.
	for bamIn in args.bam:
				
		# setup names for files
		fifo_stub = hashlib.md5(bamIn).hexdigest()
		
		# create stats object for counting hits and stuff
		bam_stats = Stats()
		
		# create hit object to track assignments
		hit_factory = HitFactory(annot=annot)
		mean_insert = mean_insert0 = 0
		
		if bamIn != "-":
			tmp = bamIn.split(".")
			fout_name = "{}.ghits".format(tmp[0])
#			if isfile(fout_name):
#				if not args.overwrite:
#					sys.stderr.write("Output file name exists ({}) and --overwrite is not specified. Skipping.".format(fout_name))
#					continue
		
		#
		# open alignments for reading. if we need to name sort them then we 
		# launch a fork to do so which will write to a fifo.
		#
		
		if args.name_sort_bam and bamIn != "-":
			# name sort the bam.  we'll fire up a sub process that calls
			# samtools sort and sort into a FIFO then open the FIFO as the
			# input to the main loop
			
			# make fifo and sort to it in a background process
			if isfile("{}.bam".format(fifo_stub)):
				os.unlink("{}.bam".format(fifo_stub))

			sys.stderr.write("+mkfifo {}.bam\n".format(fifo_stub))
			os.mkfifo("{}.bam".format(fifo_stub))

			# launch sorting process				
			sys.stderr.write("CMD: {}\n".format(cmd))
			p0 = sp.Popen(cmd.split())

			# open fifo for reading...
			p1 = sp.Popen("samtools view {}.bam".format(fifo_stub).split(), stdout=sp.PIPE)
			fin = p1.stdout

		else:
			# bam doesn't need to be sorted
			sys.stderr.write("\nNOTE: If alignments aren't name sorted then multi-mappers won't be handled.\n")

			if bamIn == "-":
				# read sam from stdin
				fin = sys.stdin
			else:
				if args.sam:
					# read sam as text
					fin = open(bamIn, "r")
				else:
					# input is bam, use samtools to read it
					p1 = sp.Popen("samtools view {}".format(bamIn).split(), stdout=sp.PIPE)
					fin = p1.stdout

		# start a timer		
		t0 = timeit.default_timer()
		
		#
		# now we loop over the alignments read in from file. bundle them into read name 
		# groups and push them out to the queue
		sys.stderr.write("parsing alignments\n")
		aln_buffer = [] # empty the buffer
		for szl in fin:
			
			if szl[0] == "@":
				# header. skip those
				continue
			
			# parse alignment
			aln = SAMAln()
			aln.parse_from_string(szl)
			num_parsed += 1
			
			if aln.aligned and aln.mapq >= args.min_mapq:
				# loop up hits for this guy
				mate = 1
				if aln.paired and aln.second_mate:
					mate = 2
				
#				total_searches += 1
#				tsearch0 = timeit.default_timer()
				aln.hits = process_bounds(aln.bounds, lktable, annot, strand_mode, aln.aln_len, args.min_overlap_ratio, mate=mate)
#				total_search_time += timeit.default_timer()-tsearch0
			
			rname = aln.qname
			if (rname != rname_last) and len(aln_buffer) > 0:
				bam_stats.parsed += 1
				if aln_buffer[0].aligned:
					bam_stats.aligned += 1
					# read is aligned, process it into final hit form
					m = process_read(aln_buffer)
					
					# update stats and update hits
					if m['result'] > 0:
						bam_stats.assigned += 1
						hit_factory.update_hits([m['hits']])
						hit_factory.update_insert_size(m['insert'])
					
					if m['result'] == ASSIGNED_GENE_UNIQUE:
						bam_stats.assigned_unique += 1
					elif m['result'] == ASSIGNED_GENE_MULTI:
						bam_stats.assigned_multi += 1
						
				if aln_buffer[0].mapq >= args.min_mapq:
					bam_stats.passed_mapq += 1
				
				# clear the buffer
				aln_buffer = []

			# append alignment				
			aln_buffer.append(aln)
			# update read name
			rname_last = rname
							
		# finished, close the file
		fin.close()
		
		if len(aln_buffer) > 0:

			bam_stats.parsed += 1
			if aln_buffer[0].aligned:
				bam_stats.aligned += 1
				# read is aligned, process it into final hit form
				m = process_read(aln_buffer)
				
				# update stats and update hits
				if m['result'] > 0:
					bam_stats.assigned += 1
					hit_factory.update_hits([m['hits']])
					hit_factory.update_insert_size(m['insert'])
				
				if m['result'] == ASSIGNED_GENE_UNIQUE:
					bam_stats.assigned_unique += 1
				elif m['result'] == ASSIGNED_GENE_MULTI:
					bam_stats.assigned_multi += 1
				
			if aln_buffer[0].mapq >= args.min_mapq:
				bam_stats.passed_mapq += 1

		# done.

		# final parsed message
		sys.stderr.write("parsed {} lines in {:0.2f} seconds\n".format(num_parsed, timeit.default_timer()-t0))

#		sys.stderr.write("\ntotal lookups: {}\ntotal lookup time: {:0.2f} seconds\n".format(total_searches, total_search_time))
		
		sys.stderr.write("\naverage fragment length: {:0.1f}\n".format(hit_factory.insert_size))

		sys.stderr.write(bam_stats.__str__())
		
		# produce output
		fout = open(fout_name, "w")
		
		fout.write("\t".join(hit_factory.get_header()) + "\n")
		
		for lout in hit_factory.get_results(min_count=0):
			fout.write("\t".join(map(str, lout)) + "\n")
		
		fout.close()
		
	return 0	

#==============================================================================
# class
#==============================================================================

#
# object to manage counters and print a summary
class Stats(object):
	def __init__(self):
		self.parsed = 0
		self.aligned = 0
		self.passed_mapq = 0
		self.assigned = 0
		self.assigned_unique = 0
		self.assigned_multi = 0
	
	def __str__(self):
		
		r1 = r2 = r3 = 0
		if self.passed_mapq > 0:
			r1 = self.assigned*100.0/self.passed_mapq
			r2 = self.assigned_unique*100.0/self.passed_mapq
			r3 = self.assigned_multi*100.0/self.passed_mapq
	
		# make a string output of the stats
		
		sout = "\nParsed {} total fragments. Of these:\n".format(self.parsed)
		sout += "  {} were aligned\n".format(self.aligned)
		sout += "  {} passed MAPQ\n".format(self.passed_mapq)
		sout += "\n"
		sout += "Of those that passed MAPQ:\n"
		sout += "  {} were assigned ({:0.2f}%)\n".format(self.assigned, r1)
		sout += "  {} were gene-unique ({:0.2f}%)\n".format(self.assigned_unique, r2)
		sout += "  {} were multi-gene ({:0.2f}%)\n".format(self.assigned_multi, r3)
		sout += "\n"
		
		return sout

#
# class to manage hits to targets and return results as a list including 
# annotation, if provided
class HitFactory(object):
	# init
	def __init__(self, annot=None):
		
		self.hits = defaultdict(float)
		self.annot = annot
		
		self.insert_size = 0
		self.insert_size_n = 0
		
		if annot is not None:
			for k in annot.keys():
				self.hits[k] = 0
	
	def update_insert_size(self, inserts):
		mu = mu0 = self.insert_size

		for s in inserts:
			self.insert_size_n += 1
			mu = mu0 + (s-mu0)*1.0/self.insert_size_n
			mu0 = mu
			
		#print "-->", self.insert_size, mu, inserts
		self.insert_size = mu
		
		return 0
			
		
	
	# input is a list of dicts where each contains a transcript_id/hit_count key/value pair
	def update_hits(self, ld):
		
		for d in ld:
			for key in d.keys():
				self.hits[key] += d[key]
		
		return 0
	
	def get_header(self):
		lout = ["transcript_id"]
		
		if self.annot is not None:
			lout += ["gene_id", "gene_name", "chrom", "strand", "length", "eff_length"]
		
		lout += ["hits", "eff_hits"]
		return lout
	
	def get_results(self, min_count=0):
		
		results = []
		eff_len = 1
		
		for tid in sorted(self.hits.keys()):
			
			if self.hits[tid] < min_count:
				continue
			
			lout = [tid]
			if self.annot is not None:
				eff_len = self.annot[tid]['length'] - self.insert_size
				lout += [ self.annot[tid]['gid'], self.annot[tid]['gname'], 
					self.annot[tid]['chrom'], self.annot[tid]['strand'], self.annot[tid]['length'], 
					eff_len ]
			
			lout.append(self.hits[tid])
			lout.append(self.hits[tid]*float(self.annot[tid]['length'])/float(eff_len))
		
			results.append(lout)
		
		return results
	
	#
	# create a summary of number of features hit and total assigned hits
	def hit_status(self):
		
		num_hit = 0
		num_hits = 0
		
		for tid in self.hits:
			if self.hits[tid] > 0:
				num_hit += 1 
				num_hits += self.hits[tid]
		
		return (num_hit, num_hits)

#
# Hit object for storing a match between a read an annotation element.
# read and target are of type 'Feature'. we can have a list of reads 
# which will occur when a read was spliced and hit the same target more than once.
class Hit(object):
	
	def __init__(self, read, target, olen):
		self.olen = olen
		self.lreads = [read]
		self.target = target
		
		self.target_name = self.target.tag
		
		self.read = None
		
		# check strand orientation. 0 for 
		# not match, 1 for match
		self.strand = 0
		if read.strand == target.strand:
			self.strand = 1
		
		self.calc_read()
	
	def __str__(self):
		
		if self.read is None:
			self.calc_read()
		
		r1 = "{}:{}-{}".format(self.read.ref, self.read.start, self.read.end)
		r2 = "{}:{}-{}".format(self.target.ref, self.target.start, self.target.end)
		sout = "{}\t{}\t{}\t{}\t{}\t{}".format(self.read.tag, self.target.tag, self.olen, self.strand, r1, r2)
		return sout
	
	#
	# create a read feature that represents all assigned read fragments
	def calc_read(self):
		
		start = end = len = 0
		
		start = self.lreads[0].start
		end = self.lreads[0].end
		
		for f in self.lreads:
			len += f.len
			start = f.start if f.start < start else start
			end = f.end if f.end > end else end
		
		self.read = Feature(ref=self.lreads[0].ref, start=start, end=end, tag=self.lreads[0].tag)
		
		return 0
	
	#
	# add a hit. this is intented for when a single read hits the same target more than once
	# as in a spliced alignment.
	def add_hit(self, h):
		# add a hit
		ftmp = Feature(ref=h.read.ref, start=h.read.start, end=h.read.end, tag=h.read.tag, strand=h.read.strand)
		# confirm that this is not identical to a region already in this hit
		keep = True
		for f in self.lreads:
			if f.is_equal(ftmp):
				keep = False
		
		if keep:
			self.lreads.append(ftmp)
			self.olen += h.olen
			self.calc_read()
			
		return 0
			
			
		

#
# look up hash of features. builds the quick-lookup hash table for the 
# annotation and looks up hits to targets from a passed Feature object.
class FeatureTable(object):
	# init
	def __init__(self, bin_size=16000):
		self.t = defaultdict(list)
		self.bin = bin_size
	
	#
	# insert a feature into the lookup table. hash the feature's position 
	# range and insert the feature into each bucket that covers its range.
	def insert(self, f):
		# hash feature
		ehash = f.hash(self.bin)
		# add to as many keys as necessary
		for key in ehash:
			self.t[key].append(f)
	
	#
	# look for overlaps of feature 'f' with any features in this lookup
	# table.
	def lookup(self, f):
		
		hits = []
		dtmp = {}

		# hash feature
		ehash = f.hash(self.bin)
		
		for key in ehash:
			if key in self.t:
				for f0 in self.t[key]:
					ovl, olen = f.overlap(f0)
					if ovl:
						h = Hit(read=f, target=f0, olen=olen)
						# plug this into a dict. we do this because a single 
						# feature may cross a bin boundary and end up hitting
						# the same target in both. this is separate from a 
						# read being spliced and having more than one 
						# feature region. that this function only deals with a 
						# single feature region.
						if f0.tag not in dtmp:
							# add the hit
							dtmp[f0.tag] = h
		
		# done
		return dtmp

# 
# feature regions are 1-based, inclusive. Represents some genomic range 
# with a strand and a tag (name). has 'overlap' function and 'is_equal'
# to determin if the feature either overlaps another or is equal to another
class Feature(object):
	# init
	def __init__(self, ref="", strand="", start=0, end=0, tag=""):
		self.ref = ref
		self.strand = strand
		self.start = start
		self.end = end
		self.len = end-start+1
		self.tag = tag
	
	def __str__(self):
		if self.tag == "":
			tag = "noname"
		else:
			tag = self.tag
			
		l = [self.ref, self.start, self.end, tag, self.strand]
		return "\t".join(map(str, l))
	
	def update(self):
		self.len = self.end-self.start+1
		
	# return length of overlap or -1 for no overlap
	def overlap(self, f):
		self.update()
		f.update()
								
		fl1 = self.len
		fl2 = f.len
		
		ovl = False
		olen = -1
		
		if (self.start <= f.end) and (self.end >= f.start):
			ovl = True
			# we have an overlap. find it's length
			l1 = self.end - f.start + 1
			l2 = f.end - self.start + 1
			
			olen = min([fl1, fl2, l1, l2])
		
		return ovl, olen
	
	#
	# return True if this and f are the same region
	def is_equal(self, f):
		
		res = False
		if self.start==f.start and self.end==f.end and self.ref==f.ref:
			res = True
		
		return res
	
	# hash feature into a bucket key which is a combination of the
	# reference name and the result of its start/end position 
	# integer divided by 'bin'
	def hash(self, bin):
		hashes = []
		start_bin = self.start/int(bin)
		end_bin = self.end/int(bin)
		
		for k in range(start_bin, end_bin+1):
			hashes.append("{}:{}".format(self.ref, k))
		
		return hashes

#
# class for a single SAM alignment. provides members for each field and
# dict access to the extended fields that come after the quality string
# the flag values are also translated into boolean members for 
# easy access
class SAMAln(object):
	# init
	def __init__(self):
		self.qname = None
		self.flag = None
		self.rname = None
		self.pos = None
		self.mapq = None
		self.cigar = None
		self.rnext = None
		self.pnext = None
		self.tlen = None
		self.seq = None
		self.qual = None
		self.ext = {}
		
		# list of type Feature containing boundaries of the alignment
		self.bounds = []
		# read length
		self.len = 0
		# aligned length (if there was soft-clipping)
		self.aln_len = 0
		
		# used to store target features this read hits from the lookup 
		# table
		self.hits = []

		self.paired = False
		self.prop_aligned = False
		self.aligned = False
		self.next_aligned = False
		self.rev = False
		self.next_rev = False
		self.first_mate = False
		self.second_mate = False
		self.secondary = False
	
	def __str__(self):
		return "{}\t{}\t{}\t{}\t{}\t{}".format(self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar)
	
	# loads alignment info from a list which would be a tab-split
	# version of the line from a sam file
	def parse_from_list(self, l):
		
		tag = val = None
		
		self.qname = l[0]
		self.flag = int(l[1])
		self.rname = l[2]
		self.pos = int(l[3])
		self.mapq = int(l[4])
		self.cigar = l[5]
		self.rnext = l[6]
		self.pnext = int(l[7])
		self.tlen = int(l[8])
		self.seq = l[9]
		self.qual = l[10]
		
		# read length
		self.len = len(self.seq)
		
		if len(l) > 11:
			for i in range(11, len(l)):
				tag, val = self.parse_optional_tag(l[i])
				self.ext[tag] = val
		
		#
		# parse the flag
		#
		
		if self.flag & 0x1:
			self.paired = True
		
		if self.flag & 0x2:
			self.prop_aligned = True
		
		if (self.flag & 0x4) == 0:
			self.aligned = True
		
		if (self.flag & 0x8) == 0:
			self.next_aligned = True
		
		if (self.flag & 0x10):
			self.rev = True
		
		if (self.flag & 0x20):
			self.next_rev = True
		
		if (self.flag & 0x40):
			self.first_mate = True
		
		if (self.flag & 0x80):
			self.second_mate = True
		
		if (self.flag & 0x100):
			self.secondary = True		
		
		#
		# if read is aligned then find the boundaries of its alignment
		if self.aligned:
			self.calc_mapped_region()
			
		# done
		return 0
	
	# 
	# parse from string just calls parse_from_list after splitting the string
	# by tabs
	def parse_from_string(self, sz):
		aln = sz.strip().split("\t")
		rres = self.parse_from_list(aln)
		return 0
	
	def is_spliced(self):
		r = re.search("[0-9]+N", self.cigar)
		if r:
			return True
		
		return False
	
	#
	# this parses the cigar and based on the alignment position it figures out 
	# the genomic start/end boundaries of the alignment. if the alignment is 
	# spliced (had a 'N' region) then it will have multiple boundaries. 
	def calc_mapped_region(self):
		# parse the cigar
		left = right = self.pos
		op_len = re.split("[A-Z]", self.cigar)
		op_typ = re.split("[0-9]+", self.cigar)
		
		# drop the extra bit that comes along with the split
		op_len = map(int, op_len[0:(len(op_len)-1)])
		op_typ = op_typ[1:len(op_typ)]
		nop = len(op_typ)
	
		# grab stuff that doesn't change
		chrom = self.rname
		strand = "-" if self.rev else "+"
	
		# loop through
		self.aln_len = 0
		for i in range(nop):
			if op_typ[i] == "M" or op_typ[i] == "D":
				right += op_len[i]
				if op_typ[i] == "M":
					self.aln_len += op_len[i]
	
			if op_typ[i] == "N":
				# start a new feature and pass the current on
				self.bounds.append(Feature(ref=chrom, strand=strand, start=left, end=right-1, tag=self.qname))
				left = right + op_len[i]
				right = left
	
		# append feature
		self.bounds.append(Feature(ref=chrom, strand=strand, start=left, end=right-1, tag=self.qname))

		return 0
	
	def parse_optional_tag(self, sz):
		tparts = sz.split(":")
		
		tag = tparts[0]
		value = tparts[2]
		if tparts[1]=="i":
			value = int(value)
		elif tparts[1] == "f":
			value = float(value)
		
		return tag, value
	
	# 
	# check if the passed SAMAln object is this alignments mate
	def is_mate(self, aln):
		
		rres = False

		if self.paired and aln.paired:		
			#if (self.rnext == aln.rname or self.rnext=="=") and (self.rname == aln.rnext):
			if (self.rnext == aln.rnext):
				if (self.pnext == aln.pos) and (self.pos == aln.pnext):
					rres = True

		return rres
	
	def has_hits(self):
		return len(self.hits) > 0
		

#==============================================================================
# defs
#==============================================================================

#
# this function handles a bundle of alignments from a single read or fragment.
# the alignments have already been intersected with the lookup table so we 
# only need to sort out what's going on with the hit.  in the case of paired
# end reads we have to pair them up first.
def process_read(item):
	
	# initial weight of alignment assigment
	wi = 1
	gid_hit = set()
	dhits = defaultdict(float)
	result = 0
	
#	annot = item[1]
#	lktable = item[2]
#	item = item[0]
	
	# need to parse the alignments into SAMAln objects. if we have paired
	# reads this also splits them into first and second mates
	aln1 = []
	aln2 = []
	pe = False
	insert_size = []
	
	for i in range(len(item)):
		aln = item[i]
				
		# skip if the read isn't aligned or if the read is paired end
		# and its mate isn't aligned
		if (not aln.aligned) or (aln.paired and (not aln.next_aligned)):
			continue
			
		if aln.paired:
			pe = True
			if aln.first_mate:
				aln1.append(aln)
			else:
				aln2.append(aln)
		
		else:
			aln1.append(aln)
	
	if pe:
		##
		# pair the alignments
		lpaired = []
		for aln in aln1:
			# find mate for this one
			for mate in aln2:
				if aln.is_mate(mate):
					lpaired.append([aln, mate])
					continue
		
		# reads are paired now...or at least those that have both mates aligned 
		# are paired.
		hits = set()
		
		if len(lpaired) > 0:
						
			# check for hits at each mate
			for i in range(len(lpaired)):
				# each alignment object has a list of the feature ranges. 
				# look up each on in the lookup table
				
				if len(lpaired[i][0].bounds)==1 and len(lpaired[i][1].bounds)==1:
					# non spliced alignment
					insert_size.append(abs(lpaired[i][0].tlen))
					#print lpaired[i][0]
					#print lpaired[i][1]
				
				# get hits from both mates
				if lpaired[i][0].has_hits() and lpaired[i][1].has_hits():
					hits1 = lpaired[i][0].hits
					hits2 = lpaired[i][1].hits
					# make a dict of the transcript ids to gene ids
					dtmp = {}
					tid1 = tid2 = []
					for m in hits1:
						tid1.append(m[0])
						dtmp[m[0]] = m[1]
					for m in hits2:
						tid2.append(m[0])
						dtmp[m[0]] = m[1]

					# now take intersection of hits1 and hits2 and update the hits set for this read
					tid_tmp = list(set(tid1).intersection(set(tid2)))
					hits.update(tid_tmp)
					# now update gid set for these transcript ids
					for tid in tid_tmp:
						gid_hit.update([dtmp[tid]])
				
				
	else:
		# single-end. look through alignments in 'aln1'
		# find hits
		hits = set()
		for aln in aln1:
			
			# for single end we'll use the read length as the insert size
			insert_size.append(aln.len)
			
			if aln.has_hits():
				for m in aln.hits:
					hits.update([m[0]])
					gid_hit.update([m[1]])
				
				
	# change hits set to a list
	hits = list(hits)
	gid_hit = list(gid_hit)
		
	# now single-end and paired-end have been reduced to a list of valid target
	# transcript ids. we can now fetch the gene ids associated with those and
	# finally assign hits
	if len(hits) > 0:
		
		result = ASSIGNED

		if len(gid_hit) > 1:
			result = ASSIGNED_GENE_MULTI
		else:
			result = ASSIGNED_GENE_UNIQUE
		
		wi = 1.0/len(hits) * 1.0/len(gid_hit)
		
		# assign hits
		for tid in hits:
			dhits[tid] += wi
		
	return {'result':result, 'hits':dhits, 'insert':insert_size }

#
# loop through all boundaries of an alignment and return targets hit
# that pass minimum overlap ratio and have the correct strand 
# orientation, if stranded

#
# strand mode: None for unstranded matching.
# FR means stranded second strand. if paired then the mates are
#    opposite one another with first mate being second strand
# RF opposite of FR.
# Return value:
# 	a dict made up of target_id/assignment_value key/value pairs.
#   a result indicating the nature of the assignment
def process_bounds(lbounds, lktable, annot, strand_mode, read_length, min_ovl_ratio, mate=1):

	keep_hits = {}
	hits = []
	for f in lbounds:
		tmp_hits = lktable.lookup(f)
		
		if not dict_empty(tmp_hits):
			
			# check the hits
			if strand_mode is not None:
				# check strand of the hits
				for k in tmp_hits.keys():
					d = tmp_hits[k]
					keep = False
					if strand_mode=="FR":
						# strands should be opposite for first mates or single end
						if mate==1 and d.strand==0:
							keep = True
						elif mate==2 and d.strand==1:
							# strand should match for second mate
							keep = True
					else:
						# same strand for first mate and single end
						if mate==1 and d.strand==1:
							keep = True
						elif mate==2 and d.strand==0:
							# second mate should be opposite strand
							keep = True

					if keep:
						if d.target_name not in keep_hits:
							keep_hits[d.target_name] = d
						else:
							keep_hits[d.target_name].add_hit(d)
			
			else:
				# unstranded - keep all of them
				
				for k in tmp_hits.keys():
					d = tmp_hits[k]
					if d.target_name not in keep_hits:
						keep_hits[d.target_name] = d
					else:
						# already in there. 
						keep_hits[d.target_name].add_hit(d)
						
	# check hits for this alignment and keep those that overlap as
	# much as the minimum ratio
	for tid in keep_hits.keys():
		if keep_hits[tid].olen*1.0/read_length > min_ovl_ratio:
			hits.append([tid, annot[tid]['gid']])
	
	return hits

#
# check if a dict is empty real quick style. returns
# false if are any contents.
def dict_empty(d):
	
	for k in d:
		return False
	
	return True
		

#
# parser of the GTF
def parse_gtf(fname):

	# variables
	gtfdb = {}
	lktable = FeatureTable(bin_size=HBIN)
	ehash = []
	szl = ""
	aln = []

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		if aln[2] != "exon":
			continue

		# parse attributes
		attr = parse_gtf_attr(aln[8])

		if attr['transcript_id'] not in gtfdb:
			if "gene_name" in attr:
				gtfdb[attr['transcript_id']] = {'gid': attr['gene_id'], 'gname': attr['gene_name'], 
					'tid': attr['transcript_id'], 'hits': 0, 'strand': aln[6], 'chrom': aln[0], 'length': 0}
			else:
				gtfdb[attr['transcript_id']] = {'gid': attr['gene_id'], 'gname': attr['gene_id'], 
					'tid': attr['transcript_id'], 'hits': 0, 'strand': aln[6], 'chrom': aln[0], 'length': 0}


		# create feature for this exon		
		ff = Feature(ref=aln[0], start=int(aln[3]), end=int(aln[4]), strand=aln[6], tag=attr['transcript_id'])
		gtfdb[attr['transcript_id']]['length'] += ff.len

		# add to hash table
		lktable.insert(ff)
	
	fin.close()

	return gtfdb, lktable

def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

#
# run a system level command. used for running alignment and samtools commands
def runcmd(cmd, verbose=True):
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	rres = os.system(cmd)
	return rres

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Quantify hits to GTF transcripts from genome alignments.")
parser.add_argument('gtf', type=str, help="GTF annotation corresponding to the reference")
# allow passing of multiple bam files at the command line
parser.add_argument('bam', type=str, metavar='bam', nargs='+', 
	help="Alignments in BAM/SAM format or - to read from stdin. If reading from stdin then input is expected to be SAM")

parser.add_argument("-p", action="store", type=int, default=1, 
	help="Number of threads for processing hits [1]")

parser.add_argument("-n", "--name-sort-bam", action="store_const", const=True, default=False,
	help="Name sort the BAM first before running quantification [off]")
parser.add_argument("-S", "--sam", action="store_const", const=True, default=False, 
	help="Input is SAM format [off]")
parser.add_argument("-l", "--min-overlap-ratio", default=0.96, type=float, action="store", 
	help="Minimum overlap as a ratio of the read, or read segment, length [0.96]")
parser.add_argument('-q', '--min-mapq', type=int, default=1, 
	help="Minimum MAPQ for counting alignment [1]")
parser.add_argument('--unweighted', action="store_const", const=True, default=False, 
	help="Do not weight hits by number of features hit or number of loci hit. Produces redundant information.")
parser.add_argument('--fr-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a reverse stranded library [off]")
parser.add_argument('--rf-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a forward stranded library [off]")
parser.add_argument('--overwrite', action="store_const", const=True, default=False, 
	help="Overwrite output file if it exists already")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

