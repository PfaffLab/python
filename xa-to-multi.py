#!/usr/bin/env python
#==============================================================================
# xa-to-multi.py
#
# Shawn Driscoll
# 20130627
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script can be used to parse BWA 'sampe' output to expand the 
# alignments to multi-alignments from the XA tags with the alignments. 
#==============================================================================

import sys, argparse, re, pysam

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

_BUFFER_SIZE = 2
_NORMAL_PAIR = 0
_SINGLETON = 2
_DISCORDANT = 1
_UNALIGNED = 3
_NOT_PAIR = -1

#==============================================================================
# main
#==============================================================================

def main(args):

	# -- variables
	paired = False
	aln_buffer = []
	rev_bit = 0x10
	cigar_sam_to_bam = { "M":0, "I":1, "D":2, "N":3, "S":4, "H":5, "P":6, "=":7, "X":8}
	
	# check input files
	if not file_exists(args.alignments):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.alignments)
		return 1


	# open the alignments as SAM or BAM
	if args.S:
		fin = pysam.Samfile(args.alignments, "r")
		fout = pysam.Samfile("dump.sam", "w", template=fin)
	else:
		fin = pysam.Samfile(args.alignments, "rb")
		fout = pysam.Samfile("dump.bam", "wb", template=fin)
	
	# loop through alignments and collect all of the junctions from alignments.
	for aln in fin:

		if aln.flag & 0x1:
			paired = True
	
		if paired:
			
			aln_buffer.append(aln)
			
			if len(aln_buffer) == _BUFFER_SIZE:
				# deal with alignments in the buffer
				pair_type = is_pair(aln_buffer)
				
				if pair_type == -1:
					# problem
					sys.stderr.write("[main] found singleton alignment (%s). abandoning the first alignment and moving on" % aln_buffer[0].qname)
					aln_buffer.pop(0)
				elif pair_type == 3:
					# unaligned pair, print back out
					for i in range(len(aln_buffer)):
						fout.write(aln_buffer[i])
					
					aln_buffer = []
				
				elif pair_type == 2:
					# we can expand on the aligned read if it has XA hits
					if aln_buffer[0].flag & 0x4:
						r1 = aln_buffer[1]
						r2 = aln_buffer[0]
					else:
						r1 = aln_buffer[0]
						r2 = aln_buffer[1]
					
					# print the original lines back out
					fout.write(aln_buffer[0])
					fout.write(aln_buffer[1])
					
					aln_buffer = []
					
					# find the XA tag
					xa_list = []
					for tag in r1.tags:
						if tag[0] == "XA":
							xa_list = tag[1].split(";")
					
					if len(xa_list) > 0:
						# deal with XA list
						print xa_list
				
				else:
					
					# print the original alignments back out
					fout.write(aln_buffer[0])
					fout.write(aln_buffer[1])
					
					
					# find the XA tags
					xa_list_1 = []
					xa_list_2 = []

					# ENSMUST00000163854,+3517,101M,1;
					
					temp = aln_buffer[0].opt("XA")
					if temp:
						xa_list = temp.split(";")
						# make a list of alignments from the XA alignments
						alns_1 = []
						
						for alt in xa_list:
							if len(alt) > 0:
								ltmp = alt.split(",")
								ptmp = list(ltmp[1])
								
								# build cigar from string
								cigar = []						
								ops, lengths = parse_cigar(ltmp[2])
								for i in range(len(ops)):
									cigar.append((cigar_sam_to_bam[ops[i]], lengths[i]))
								
								# create alignment object
								b = pysam.AlignedRead()
								b.qname = aln_buffer[0].qname
								b.seq = aln_buffer[0].seq
								b.flag = aln_buffer[0].flag
								b.flag |= 0x100
								
								if ptmp[0] == "-":
									b.flag |= rev_bit
								else:
									b.flag &= ~rev_bit
								
								b.tid = fout.gettid(ltmp[0])
								
								b.pos = int("".join(ptmp[1:]))
								b.mapq = 255
								b.cigar = cigar
								b.mrnm = 0
								b.mpos = 0
								b.isize = 0
								b.qual = aln_buffer[0].qual
								b.tags = [ ("NM", int(ltmp[-1])) ]
								alns_1.append(b)
							

					temp = aln_buffer[1].opt("XA")
					if temp:
						xa_list = temp.split(";")
						# make a list of alignments from the XA alignments
						alns_2 = []
						
						for alt in xa_list:
							if len(alt) > 0:
								ltmp = alt.split(",")
								ptmp = list(ltmp[1])
								
								# build cigar from string
								cigar = []						
								ops, lengths = parse_cigar(ltmp[2])
								for i in range(len(ops)):
									cigar.append((cigar_sam_to_bam[ops[i]], lengths[i]))
						
								
								# create alignment object
								b = pysam.AlignedRead()
								b.qname = aln_buffer[1].qname
								b.seq = aln_buffer[1].seq
								b.flag = aln_buffer[1].flag
								b.flag |= 0x100
								
								if ptmp[0] == "-":
									b.flag |= rev_bit
								else:
									b.flag &= ~rev_bit
								
								b.tid = fout.gettid(ltmp[0])
								
								b.pos = int("".join(ptmp[1:]))
								b.mapq = 255
								b.cigar = cigar
								b.mrnm = 0
								b.mpos = 0
								b.isize = 0
								b.qual = aln_buffer[1].qual
								b.tags = [ ("NM", int(ltmp[-1])) ]
								alns_2.append(b)
						
					
					print alns_1[0]
					print alns_2[0]
					
					print alns_2[0].aend

					aln_buffer = []
		
		else:
			if aln.flag & 0x4:
				# print unaligned back out
				fout.write(aln)
			
			else:
				
				# deal with this alignment
				print aln.tags 
					
		
	
	fin.close()
	fout.close()

	return 0

def is_pair(laln):
	
	if laln[0].qname != laln[1].qname:
		# not a pair...different query names
		return _NOT_PAIR

	if (laln[0].flag & 0x4) and (laln[1].flag & 0x4):
		# unaligned pair
		return _UNALIGNED

	if ((laln[0].flag & 0x100) != 0) and ((laln[1].flag & 0x100) == 0):
		sys.stderr.write("[is_pair] warning: secondary and primary alignments mixed together\n")
		return _NOT_PAIR

	if ((laln[0].flag & 0x100) == 0) and ((laln[1].flag & 0x100) != 0):
		sys.stderr.write("[is_pair] warning: secondary and primary alignments mixed together\n")
		return _NOT_PAIR
	
	if (laln[0].flag & 0x8) or (laln[1].flag & 0x8):
		# singleton
		return _SINGLETON
	
	if laln[0].tid != laln[1].tid:
		# discordant pair
		return _DISCORDANT
	
	if (laln[0].pos != laln[1].pnext) or (laln[1].pos != laln[0].pnext):
		sys.stderr.write("[is_pair] Warning: possible messed up mate fields\n")
	
	# return 0 for normal pairs
	return _NORMAL_PAIR

def parse_cigar(sz):
	#
	# parse CIGAR notation into two lists. one list of CIGAR types and another
	# of the lengths associated with the type
	#
	split_1 = re.split("[A-Z]",sz)
	split_2 = re.split("[0-9]+",sz)

	split_1 = split_1[0:(len(split_1)-1)]
	split_2 = split_2[1:]

	return (split_2,map(int,split_1))

def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Builds a junction library from BAM/SAM alignments.")
parser.add_argument('alignments', type=str, help="BAM alignments (or SAM with -S)")
parser.add_argument("-S", dest="S", action="store_const", const=True, default=False, 
				help="Alignments are SAM (default expected: BAM)")
#parser.add_argument("-o", dest="o", action="store", default="junctions", 
#				help="Stub for output files (default: junctions)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
