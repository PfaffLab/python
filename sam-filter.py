#!/usr/bin/env python
#==============================================================================
# sam-filter.py
#
# Shawn Driscoll
# 20120510
#
# Filter SAM data based on certain flags (typical of bowtie2/tophat output)
#==============================================================================

import sys,argparse,re

#==============================================================================
# parse args
#==============================================================================

parser = argparse.ArgumentParser(description="Filter SAM file")
parser.add_argument("infile",action="store",type=str,help="FASTQ input or - to read from stdin")
parser.add_argument("-xm",action="store",dest="xmf",type=int,default=-1,help="Filter by XM flag. Discard alignments with value greater than this (default: off)")
parser.add_argument("-asl",action="store",dest="asl",type=int,default=-1,help="Filter by alignment score. Discard alignments with value less than this (default: off)")
parser.add_argument("-mq",action="store",dest="mq",type=int,default=-1,help="Filter by MAPQ. Discard reads with MAPQ lower than this. (default: off)")
args = parser.parse_args()

#==============================================================================
# main scirpt
#==============================================================================

if args.infile == "-":
	in1 = sys.stdin
else:
	in1 = open(args.infile,"r")

nReads = 0
nFiltered = 0
bKeep = True

for read in in1:
	read = read.strip()
	
	nReads += 1
	bKeep = True

	if args.xmf >= 0:
		m = re.search('XM\:i\:([0-9]+)',read)
		if m:
			if int(m.group(1)) > args.xmf:
				bKeep = False

	if args.asl >= 0:
		m = re.search('AS\:i\:([\-0-9]+)',read)
		if m:
			if int(m.group(1)) < args.asl:
				bKeep = False

	if args.mq > 0:
		arl = read.split("\t")
		if int(arl[4]) < args.mq:
			bKeep = False

	if bKeep:
		sys.stdout.write(read + "\n")
	else:
		nFiltered += 1


sys.stderr.write("> filtered " + str(nFiltered) + " reads\n")

