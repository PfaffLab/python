#!/usr/bin/env python
#==============================================================================
# bam-genome-bin-counts.py
#
# Shawn Driscoll
# 20160219
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Counts alignments to a genome into bins of user specified size
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib
from multiprocessing import cpu_count

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

	read_num_align = 0 # used for counting number of alignments per read

	read_strand_align = []

	current_read = ""
	last_read = ""
	paired = False
	base_weight = 1
	hits = []
	annot = {}
	lktable = {}
	fifo_stub = hashlib.md5(args.bam).hexdigest()
	max_threads = cpu_count()/2


	# -------------------------------------------------------------------------
	# load GTF annotation
	# -------------------------------------------------------------------------

	sys.stderr.write("loading FAI and building lookup table for genome...\n")
	dlk_table, lbin_id = parse_fai(args.fai, args.bin_size)
	sys.stderr.write("created {} bins\n".format(len(lbin_id)))


	# -------------------------------------------------------------------------
	# ready to parse the bam
	# -------------------------------------------------------------------------

	if args.name_sort_bam:
		# name sort the bam.  we'll fire up a sub process that calls
		# samtools sort and sort into a FIFO then open the FIFO as the
		# input to the main loop
		
		# make fifo and sort to it in a background process
		if isfile("{}.bam".format(fifo_stub)):
			os.unlink("{}.bam".format(fifo_stub))

		sys.stderr.write("+mkfifo {}.bam\n".format(fifo_stub))
		os.mkfifo("{}.bam".format(fifo_stub))

		# launch sorting process
		cmd = "samtools sort -n -@ {} {} {}".format(max_threads, args.bam, fifo_stub)
		sys.stderr.write("CMD: {}\n".format(cmd))
		p0 = sp.Popen(cmd.split())

		# open fifo for reading...
		p1 = sp.Popen("samtools view {}.bam".format(fifo_stub).split(), stdout=sp.PIPE)
		fin = p1.stdout


	else:
		# bam doesn't need to be sorted
		sys.stderr.write("I hope your alignments are already name sorted...\n")

		if args.bam == "-":
			# read sam from stdin
			fin = sys.stdin
		else:
			if args.sam:
				# read sam as text
				fin = open(args.bam, "r")
			else:
				# input is bam, use samtools to read it
				p1 = sp.Popen("samtools view {}".format(args.bam).split(), stdout=sp.PIPE)
				fin = p1.stdout

	# -- MAIN LOOP

	sys.stderr.write("starting main loop...\n")

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

				# count the number of references the read aligns to
				ref_set = set()
				for k in hits:
					ref_set.update([dlk_table[k]['chrom']])

				weight = base_weight/read_num_align * 1.0/(len(ref_set)**2)

				# loop through hits
				for i in range(len(hits)):
					bid = hits[i]
					if read_strand_align[i]==STRAND_POS:
						dlk_table[bid]['phits'] += weight
					else:
						dlk_table[bid]['nhits'] += weight

			# reset read bins
			hits = []
			read_num_align = 0
			read_strand_align = []


		# 
		# check alignment to see if we want to count it
		#

		if not (aln[1] & 0x4 or aln[1] & 0x8):
			num_aligns += 1
			aln_ok = True

			if int(aln[4]) < args.min_mapq:
				aln_ok = False
			elif paired and aln[6] != "*":
				aln_ok = False

			# check this dog out
			if aln_ok:
				num_aligns_accepted += 1
				read_num_align += 1

				aln_feature = aln_to_feature(aln)

				# everything is a hit
				for f in aln_feature:
					hits += hash_feature(f, args.bin_size)

					for i in range(len(hits)):
						if paired:
							if (aln[1] & 0x40):
								if (aln[1] & 0x10):
									# reverse strand
									read_strand_align.append(STRAND_NEG)
								else:
									read_strand_align.append(STRAND_POS)
							else:
								# second mate
								if (aln[1] & 0x10):
									# reverse strand
									read_strand_align.append(STRAND_POS)
								else:
									read_strand_align.append(STRAND_NEG)
						else:
							# not paired
							if (aln[1] & 0x10):
								read_strand_align.append(STRAND_NEG)
							else:
								read_strand_align.append(STRAND_POS)

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
		pass
	
	num_reads += 1
	if len(hits) > 0:
		num_reads_assigned += 1
		# get ready to assign hits

		# count the number of references the read aligns to
		ref_set = set()
		for k in hits:
			ref_set.update([dlk_table[k]['chrom']])

		weight = base_weight/read_num_align * 1.0/(len(ref_set)**2)

		# loop through hits
		for i in range(len(hits)):
			bid = hits[i]
			if read_strand_align[i]==STRAND_POS:
				dlk_table[bid]['phits'] += weight
			else:
				dlk_table[bid]['nhits'] += weight

	sys.stderr.write("parsed/reads/assigned = {}/{}/{:0.1f}%\n".format(num_parsed, 
		num_reads, num_reads_assigned*100.0/num_reads))

	print_report(num_parsed, num_aligns, num_aligns_accepted, num_unaligned, num_reads, num_reads_assigned)

	# produce results
	# print header...
	print "\t".join(["chrom", "start", "end", "bin_id", "pstrand_hits", "nstrand_hits", "hits"])
	for i in range(len(lbin_id)):
		bid = lbin_id[i][0]
		db = dlk_table[bid]

		lout = [db['chrom'], db['start'], db['end'], db['id'], db['phits'], db['nhits'], db['phits']+db['nhits']]

		print "\t".join(map(str, lout))


	return 0

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
# parse_fai
# parse the genome FAI file and make all of the buckets
def parse_fai(fai, bin_size):

	dbins = {}
	lbin_id = []
	arl = []
	bin_id = ""
	bin_idx = 0
	fin = open(fai, "r")

	for szl in fin:
		arl = szl.strip().split("\t")
		refname = arl[0]
		reflen = int(arl[1])
		num_bins = reflen/bin_size
		for i in range(num_bins+1):
			bin_hash = "{}:{}".format(refname, i)
			bin_id = "GBIN_{:08d}".format(bin_idx)
			lbin_id.append([bin_hash, bin_id])
			bin_idx += 1
			dbins[bin_hash] = {"id": bin_id, "chrom":refname, "start": i*bin_size+1, "end": i*bin_size+bin_size, "phits": 0, "nhits": 0 }

	fin.close()
	return dbins, lbin_id


# --
# check_lktable
# used for looking up hits to features in the lookup table. returns
# a list of transcript ids hit by feature
def check_lktable(feature, lktable):

	fhahes = []
	hits = []

	# hash the feature
	fhashes = hash_feature(feature)

	for k in fhashes:
		hits.append(k)

	return hits

# --
# hash_feature
# feature is a list with three elements: chrom, start, end
def hash_feature(lfeature, bin_size):

	hashes = []
	start_bin = lfeature[1]/bin_size
	end_bin = lfeature[2]/bin_size

	for k in range(start_bin, end_bin+1):
		hashes.append("{}:{}".format(lfeature[0], k))

	return hashes


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

parser = argparse.ArgumentParser(description="Quantfies reads into genome bins.")
parser.add_argument('fai', type=str, help="FAI index for the genome FASTA the alignments are relative to")

parser.add_argument('bam', type=str, 
	help="Alignments in BAM/SAM format or - to read from stdin. If reading from stdin then input is expected to be SAM")

parser.add_argument("-n", "--name-sort-bam", action="store_const", const=True, default=False,
	help="Name sort the BAM first before running quantification [off]")

parser.add_argument("-S", "--sam", action="store_const", const=True, default=False, 
	help="Input is SAM format [off]")

parser.add_argument('-q', '--min-mapq', type=int, default=1, 
	help="Minimum MAPQ for counting alignment [1]")

parser.add_argument("-b", "--bin-size", type=int, default=16000, action="store", 
	help="Bin size for alignment counting [16000]")

parser.add_argument("-o", "--overlapping-bins", action="store_const", const=True, default=False, 
	help="If set then bins will be overlapping bins of 1/2 the length. This doubles the total number of bins and may provide smoother information [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

