#!/usr/bin/env python
#==============================================================================
# bam-gcounts.py
#
# Shawn Driscoll
# 20160219
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This is one of those things that counts hits to features annotated in a GTF
# from read alignments to a genome.
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib
import multiprocessing as mp

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
HBIN = 16000
STRAND_POS = "+"
STRAND_NEG = "-"

#==============================================================================
# main
#==============================================================================

def main(args):
	# variables
	total_depth = 0
	nlines = 0
	
	# variables for counting things
	num_parsed = 0 # this should be total lines from the bam
	num_aligns = 0 # this will be total lines in the bam that were aligned
	num_aligns_accepted = 0 # total alignments that pass the filters to be counted
	num_unaligned = 0 # total unaligned. total distinct reads aligned is num_reads-num_unaligned
	num_reads = 0 # total distinct reads
	num_reads_assigned = 0 # total distinct reads assigned

	current_read = ""
	last_read = ""
	paired = False
	base_weight = 1
	hits = []
	annot = {}
	lktable = {}
	fifo_stub = ""
	max_threads = cpu_count()/2

	bamIn = ""
	foutName = ""
	tmp = []


	# -------------------------------------------------------------------------
	# load GTF annotation
	# -------------------------------------------------------------------------

	sys.stderr.write("loading GTF annotation...\n")
	annot, lktable = parse_gtf(args.gtf)
	sys.stderr.write("found {} features\n".format(len(annot.keys())))


	# -------------------------------------------------------------------------
	# ready to parse the bam
	# -------------------------------------------------------------------------

	for bamIn in args.bam:
		sys.stderr.write("\nProcessing {}\n".format(bamIn))

		# reset counting variables
		num_parsed = 0 
		num_aligns = 0 
		num_aligns_accepted = 0 
		num_unaligned = 0 
		num_reads = 0 
		num_reads_assigned = 0 
		foutName = None

		# reset count hash
		reset_count_table(annot)

		#
		# make fifo name
		fifo_stub = hashlib.md5(bamIn).hexdigest()

		#
		# create output file name
		if bamIn != "-":
			tmp = bamIn.split(".")
			foutName = "{}.ghits".format(tmp[0])
			if isfile(foutName):
				sys.stderr.write("Output file name exists ({}). Skipping".format(foutName))
				continue

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
			cmd = "samtools sort -n -@ {} -o {} {}".format(max_threads, "{}.bam".format(fifo_stub), bamIn)
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

		# -- MAIN LOOP

		for szl in fin:
			# skip header lines
			if szl[0] == "@":
				continue

			# count parsed lines
			num_parsed += 1

			aln = szl.strip().split("\t")
			aln[1] = int(aln[1])
			base_weight = 1
			if aln[1] & 0x1:
				base_weight = 0.5
				paired = True

			# get current read name
			current_read = aln[0]

			# does this read match the last?
			if current_read != last_read:
				num_reads += 1
				if len(hits) > 0:
					num_reads_assigned += 1
					# get ready to assign hits

					# get gene id set
					gset = tid_to_gset(hits, annot)
					weight = float(base_weight)/len(hits) * 1.0/(len(gset)**2)

					# loop through hits
					for tid in hits:
						annot[tid]['hits'] += weight

				hits = []

			if not (aln[1] & 0x4 or aln[1] & 0x8):
				num_aligns += 1
				# check this dog out
				if int(aln[4]) >= args.min_mapq:
					num_aligns_accepted += 1
					aln_feature = aln_to_feature(aln)
					hits0 = []
					for f in aln_feature:
						hits0 += check_lktable(f, lktable)

					if args.r_stranded:
						hits += check_hit_rstrand(aln, hits0, annot)
					elif args.f_stranded:
						hits += check_hit_fstrand(aln, hits0, annot)
					else:
						hits += hits0

			else:
				num_unaligned += 1

			last_read = current_read

			if (num_parsed % 1000000) == 0:
				sys.stderr.write("parsed/reads/assigned = {}/{}/{:0.1f}%\n".format(num_parsed, 
					num_reads, num_reads_assigned*100.0/num_reads))

		fin.close()
		# remove the fifo
		try:
			os.unlink("{}.bam".format(fifo_stub))
		except:
			sys.stderr.write("no temp sorting file found (so I won't delete it)\n")
		
		num_reads += 1
		if len(hits) > 0:
			# get ready to assign hits
			num_reads_assigned += 1
			# get gene id set
			gset = tid_to_gset(hits, annot)
			weight = float(base_weight)/len(set(hits)) * 1.0/(len(gset)**2)

			# loop through hits
			for tid in hits:
				annot[tid]['hits'] += weight

		sys.stderr.write("\nparsed/reads/assigned = {}/{}/{:0.1f}%\n".format(num_parsed, 
			num_reads, num_reads_assigned*100.0/num_reads))

		print_report(num_parsed, num_aligns, num_aligns_accepted, num_unaligned, num_reads, num_reads_assigned)

		# produce results
		# print header...
		if foutName is not None:
			fout = open(foutName, "w")
		else:
			fout = sys.stdout
			
		fout.write("\t".join(["transcript_id", "gene_id", "gene_name", "chrom", "strand", "length", "count"]))
		fout.write("\n")
		for tid in sorted(annot.keys()):
			lout = [tid, annot[tid]['gid'], annot[tid]['gname'], 
				annot[tid]['chrom'], annot[tid]['strand'], annot[tid]['length'], 
				annot[tid]['hits']]
			fout.write("\t".join(map(str, lout)))
			fout.write("\n")

		if foutName is not None:
			fout.close()


	return 0


class WorkQueue(object):
	#
	# initialization
	def __init__(self, threads=1):
		self.threads = threads
		
		# the comparison of alignments to the hashed index will happen 
		# in this class so this class needs have the lookup table
		self.lktable = None
		# annotation to relate 
		self.annot = None
	
		# the queue
		self.q = Queue()
		
		# we will end up with one of these from each thread. we'll merge them 
		# down later
		self.bins = {}
	
		#
		# create the threads
		for i in range(self.threads):
			t = threading.Thread(target=self.worker)
			t.daemon = True
			t.start()
						
	#
	# main worker. this receives the input from the queue and handles it. 
	# each item passed to a worker will be a list of alignments all from
	# the same read.	
	def worker(self):
		name = threading.currentThread().getName()
		
		# create a dict to hold the results from this thread
		if name not in self.bins:
			self.bins[name] = {}

		while True:
			item = self.q.get()
			
			# 'item' will be a list of alignments for a single read name. 
			# look up hits and deal with it then assign to those features
			# in this thread's dict in 'self.bins'
			
			self.q.task_done()		

	#
	# this function merges the separate thread bins down into a single 
	# dict of counts
	def merge_bins(self):
		
		tnames = self.bins.keys()
		k0 = tnames[0]
		tid = ""
		
		for k in tnames[1:len(tnames)]:
			for tid in self.bins[k].keys():
				# add entry if necessary
				if tid not in self.bins[k0]:
					self.bins[k0][tid] = 0
				
				# increment count
				self.bins[k0][tid] += self.bins[k][tid]
		
		# done, return the merged bin
		return self.bins[k0]
		
		
		
		
# --
# print_report
#
def print_report(np, na, naa, nu, nr, nra):
	sys.stderr.write("""
total rows parsed:           {}
total alignments:            {}
total accepted alignments:   {} ({:0.1f}% of alignments)
total distinct reads:        {}
total reads aligned:         {} ({:0.1f}% of reads)
total reads assigned:        {} ({:0.1f}% of aligned)
average alignments per read: {:0.1f}

""".format(np, na, naa, naa*100.0/na, nr, nr-nu, (nr-nu)*100.0/nr, 
		nra, nra*100.0/(nr-nu), na*1.0/(nr-nu)))

	return 0

# --
# check_hit_rstrand
# if alignments are from a stranded library then we have to do this.
def check_hit_rstrand(aln, hits, annot):

	# update list with valid hits based on strand
	hits0 = []
	# get strand of alignment
	aln_strand = STRAND_POS
	if aln[1] & 0x10:
		aln_strand = STRAND_NEG

	for tid in hits:
		tstrand = annot[tid]['strand']

		if aln[1] & 0x1:
			if aln[1] & 0x40:
				# paird and first mate - strand should be opposite
				if aln_strand != tstrand:
					hits0.append(tid)

			elif aln[1] & 0x80:
				if aln_strand == tstrand:
					hits0.append(tid)
		else:
		 	# single end
		 	if aln_strand != tstrand:
		 		hits0.append(tid)

	return hits0

# --
# check_hit_fstrand
# if alignments are from a stranded library then we have to do this.
def check_hit_fstrand(aln, hits, annot):

	# update list with valid hits based on strand
	hits0 = []
	# get strand of alignment
	aln_strand = STRAND_POS
	if aln[1] & 0x10:
		aln_strand = STRAND_NEG

	for tid in hits:
		tstrand = annot[tid]['strand']

		if aln[1] & 0x1:
			if aln[1] & 0x40:
				# paird and first mate - strand should be opposite
				if aln_strand == tstrand:
					hits0.append(tid)
			elif aln[1] & 0x80:
				if aln_strand != tstrand:
					hits0.append(tid)
		else:
		 	# single end
		 	if aln_strand == tstrand:
		 		hits0.append(tid)

	return hits0


#
# parser of the GTF
def parse_gtf(fname):
	# variables
	gtfdb = {}
	lktable = {}
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

		
		# increment length
		gtfdb[attr['transcript_id']]['length'] += (int(aln[4])-int(aln[3]))+1

		# hash the feature
		lfeature = [aln[0], int(aln[3]), int(aln[4]), attr['transcript_id'], aln[6]]
		ehash = hash_feature(lfeature)
		for k in ehash:
			if k not in lktable:
				lktable[k] = []
			lktable[k].append(list(lfeature))

	fin.close()

	return gtfdb, lktable

#
# function to reset the count table between files
def reset_count_table(ll):
	for tid in ll.keys():
		ll[tid]['hits'] = 0
	
	return 0


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

# --
# check_lktable
# used for looking up hits to features in the lookup table. returns
# a list of transcript ids hit by feature
def check_lktable(feature, lktable):

	fhahes = []
	hits = []
	# length of the feature
	flen = feature[2]-feature[1]+1

	# hash the feature
	fhashes = hash_feature(feature)

	for k in fhashes:
		if k in lktable:
			for f in lktable[k]:
				# check hit to f confirming the overlap is at least
				# xx% of the feature length
				if feature_overlap(feature, f) > (flen*0.96):
					# append transcript id
					hits.append(f[3])

	return hits

# --
# hash_feature
# feature is a list with three elements: chrom, start, end
def hash_feature(lfeature):

	hashes = []
	start_bin = lfeature[1]/HBIN
	end_bin = lfeature[2]/HBIN

	for k in range(start_bin, end_bin+1):
		hashes.append("{}:{}".format(lfeature[0], k))

	return hashes

# --
# tid_to_gset
# look up the gene ids associated with the passed list
# of transcript ids
def tid_to_gset(ltid, annot):
	gset = set()

	for tid in ltid:
		gset.update([annot[tid]['gid']])

	return(list(gset))

# --
# feature_overlap
# returns length of the overlap of the two features.  expects the same list
# input for a feature as hash_feature [chrom, start, end]
def feature_overlap(f1, f2):
	rres = 0
	
	f1_len = f1[2]-f1[1]+1
	f2_len = f2[2]-f2[1]+1
	
	if f1[1] <= f2[2] and f1[2] >= f2[1]:

		# [============]
		#    [======]

		# what is the overlap length?
		l1 = f1[2]-f2[1]+1
		l2 = f2[2]-f1[1]+1

		# start with min of the above two lengths
		rres = min([l1, l2, f1_len, f2_len])


	return rres

# --
# aln_to_feature
# converts sam alignment into feature regions.
def aln_to_feature(aln):
	# vars

	cigar = aln[5]
	op_len = re.split("[A-Z]", cigar)
	op_typ = re.split("[0-9]+", cigar)
	lfeat = []

	# drop the extra bit that comes along with the split
	op_len = map(int, op_len[0:(len(op_len)-1)])
	op_typ = op_typ[1:len(op_typ)]
	nop = len(op_typ)

	# grab stuff that doesn't change
	chrom = aln[2]
	strand = STRAND_POS
	if int(aln[1]) & 0x10:
		strand = STRAND_NEG


	# loop through

	left = int(aln[3])
	right = left
	for i in range(nop):
		if op_typ[i] == "M" or op_typ[i] == "D":
			right += op_len[i]

		if op_typ[i] == "N":
			# start a new feature and pass the current on
			lfeat.append([chrom, left, right, strand])
			left = right + op_len[i]
			right = left

	# append feature
	lfeat.append([chrom, left, right, strand])

	return lfeat



def runcmd(cmd, returnProcess):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Quantfies reads against a reference using kmers.")
parser.add_argument('gtf', type=str, help="GTF annotation corresponding to the reference")

# allow passing of multiple bam files at the command line
parser.add_argument('bam', type=str, metavar='bam', nargs='+', 
	help="Alignments in BAM/SAM format or - to read from stdin. If reading from stdin then input is expected to be SAM")
parser.add_argument("-n", "--name-sort-bam", action="store_const", const=True, default=False,
	help="Name sort the BAM first before running quantification [off]")
parser.add_argument("-S", "--sam", action="store_const", const=True, default=False, 
	help="Input is SAM format [off]")
parser.add_argument("-l", "--min-overlap-ratio", default=0.96, type=float, action="store", 
	help="Minimum overlap as a ratio of the read, or read segment, length [0.96]")
parser.add_argument('-q', '--min-mapq', type=int, default=1, 
	help="Minimum MAPQ for counting alignment [1]")
parser.add_argument('--r-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a reverse stranded library [off]")
parser.add_argument('--f-stranded', action="store_const", const=True, default=False, 
	help="Alignments are from a forward stranded library [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

