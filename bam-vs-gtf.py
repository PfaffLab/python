#!/usr/bin/env python
#==============================================================================
# bam-vs-gtf.py
#
# Shawn Driscoll
# 20130625
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Intersects alignments with a GTF (a simple read counter)
#==============================================================================

import sys, argparse, pysam

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

_CHROM = 0
_START = 1
_END = 2
_STRAND = 3
_TYPE = 4
_ATTR = 5

HASH_BIN = 12000

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	szl = ""
	lk_table = {}
	paired = False
	hit_ratio = 1
	
	counted = 0
	skipped = 0
	unaligned = 0
	ambiguous = 0
	multi_hit = 0
	multi_aligned = 0
	ppos = 0
	

	#----------- check input files

	if not file_exists(args.alignments):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.alignments)
		return 1
	if not file_exists(args.gtf):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.gtf)
		return 1

	#----------- get to work

	sys.stderr.write("[main] parsing GTF file...\n")
	lk_table = gtf_lookup_hash(args.gtf_file)

	# open snp file and find stuff

	sys.stderr.write("[main] intersecting SNPs with GTF features...\n")

	if args.S:
		sf = pysam.Samfile(args.alignments,"r")
	else:
		sf = pysam.Samfile(args.alignments, "rb")

		
	if args.u:
		# unsorted alignments so we don't need to use a buffer
		
		for aln in sf:
			
			paired = aln.flag & 0x1
			hit_ratio = 1
						
			if alignment_accepted(aln, args):
				
				# loop through cigar ops and intersect with gtf
				hits = []
				ppos = aln.pos
				for op in aln.cigar:
					
					if op[0] == 0:
						# make a list of the target reference and the coordinates of this segment of the
						# alignment
						lpos = [sf.header['SQ'][aln.tid]['SN'], ppos+1, ppos+op[1]]
						hits += lookup_position(lk_table, args, lpos)
						ppos += op[1]
					elif op[0] == 2 or op[0] == 3:
						# - N or D, advance the position but don't count as overlap
						ppos += op[1]
					
				if len(hits) > 0:
					
					# adjust hit_ratio for singletons, if paired and if allowed
					if paired and (aln.flag & 0x8) and not args.no_singletons:
						hit_ratio = 1
					elif paired:
						hit_ratio = 0.5
					
					hit_ratio /= float(len(hits))
						
				
			else:
				skipped += 1
			
			
		
	else:
		# sorted by name so we'll use a buffer to make sure we catch multi-alignments
		# properly
		pass


	return 0


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

# 
# parse_gtf_attr
def parse_gtf_attr(field):

	attrs = {}

	# split into fields
	l1 = field.split(";")

	# split each field into key and value
	l2 = []
	for i in range(len(l1)):
		l2.append(l1[i].split("\""))


	for i in range(len(l2)):
		if len(l2[i]) > 1:
			attrs[l2[i][0].strip()] = l2[i][1].strip()

	return(attrs)		

#
# make a fast lookup hash of the GTF
#
def gtf_lookup_hash(fname):

	lk = {}

	#- open file -#

	fin = open(fname, "r")

	#- parse file and build hash

	for szl in fin:
		ll = szl.strip().split("\t")
		
		if ll[2] == "exon":
			kid = hash_pos(ll[0], ll[3])

			if kid not in lk:
				lk[kid] = []
						
			lk[kid].append([ll[0], int(ll[3]), int(ll[4]), ll[6], ll[2], parse_gtf_attr(ll[8])])

	fin.close()

	return(lk)

def hash_pos(rname, pos):
	mod = 16000
	bin = int(pos)/mod
	kid = rname + str(bin)
	return kid

# 
# look up a position in the hash, return list of feature names that intersect
#
def lookup_position(lk, args, lpos):

	result = []
	llen = int(lpos[2])-int(lpos[1])+1
	
	# hash both start and end
	query_id1 = hash_pos(lpos[0], lpos[1])
	query_id2 = hash_pos(lpos[0], lpos[2]) 

	if query_id1 in lk:
		# check for intersections with features in this bucket
		ll = lk[query_id1]
		
		for i in range(len(ll)):
			lf = ll[i]
			
			# if int(lpos[1]) <= int(lf[_END]) and int(lpos[2]) >= int(lf[_START]):
			min_overlap = int(lf[_END])-int(lpos[1])+1
			min_overlap = min(min_overlap, int(lpos[2])-int(lf[_START])+1)
			
			if min_overlap >= args.m:
				# got it
				result.append(lf[_ATTR][args.i])
	
	elif query_id2 in lk or query_id2 != query_id1:
		# check for intersections with features in this bucket
		ll = lk[query_id2]
		
		for i in range(len(ll)):
			lf = ll[i]
			
			# if int(lpos[1]) <= int(lf[_END]) and int(lpos[2]) >= int(lf[_START]):
			min_overlap = int(lf[_END])-int(lpos[1])+1
			min_overlap = min(min_overlap, int(lpos[2])-int(lf[_START])+1)
			
			if min_overlap >= args.m:
				# got it
				result.append(lf[_ATTR][args.i])
		
	return result

#
# alignment accepted
#
def alignment_accepted(aln, args):
	
	result = True
	
	
	
	return result

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Simple read counter for SAM/BAM alignments verses a GTF annotation")
parser.add_argument('alignments', type=str, help="BAM or SAM alignments (use -S if SAM)")
parser.add_argument('gtf', type=str, help="GTF reference file")
parser.add_argumnet('-i', type=str, default='gene_id', dest="i", help="GTF attribute for summing hits (default: gene_id)")
parser.add_argument('-m', type=int, default=1, dest="m", help="Minimum overlap to count as a hit (default: 1)")
parser.add_argument('-n', action="store_const", const=True, default=False, dest="n", 
				help="Discard ambiguous hits (default: count at each)")
parser.add_argument('-S', dest="S", action="store_const", const=True, default=False, help="Input is SAM (default: bam)")
parser.add_argument('-q', action="store", dest="q", default=0, help="Minimum MAPQ to quantify (default: 0)")
parser.add_argument('--count-secondary', action="store_const", dest="count_secondary", const=True, default=False,
				help="Count secondary alignments (by the 0x100 flag) as hits (default: off)")
parser.add_argument('--no-discordant', dest="no_discordant", action="store_const", const=True, default=False, 
				help="Do not quantify discordant paired end alignments (default: off)")
parser.add_argument('--no-singletons', dest="no_singletons", action="store_const", const=True, default=False,
				help="Do not quantify singleton paired-end alignments (default: off)")
parser.add_argument('-u', action="store_const", dest="u", const=True, default=False, 
				help="The input alignments are not name sorted. This means if a read appears more than once in the alignments and is not flagged as secondary it will be counted more than once.")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
