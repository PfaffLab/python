#!/usr/bin/env python
#==============================================================================
# fastq-subset.py
#
# Shawn Driscoll
# 20120502
#
# base loop code provided by Simon Anders.
# writes output to outfile.
#==============================================================================

import sys,random,itertools
import argparse
import HTSeq

#==============================================================================
# parse args
#==============================================================================

parser = argparse.ArgumentParser(description="Select randomly select reads from FASTQ data by percentage of total reads in file")
parser.add_argument("fraction",type=float,action="store",help="Fraction of reads as a decimal (range: 0,1)")
parser.add_argument("outfile",action="store",type=str,help="name for FASTQ output file")
parser.add_argument("infile",action="store",type=str,help="FASTQ input or - to read from stdin")
args = parser.parse_args()

#==============================================================================
# main scirpt
#==============================================================================

if args.infile == "-":
	in1 = iter(HTSeq.FastqReader(sys.stdin))
else:
	in1 = iter(HTSeq.FastqReader(args.infile))

out1 = open(args.outfile, "w")

for read in in1:
	if random.random() < args.fraction:
		read.write_to_fastq_file(out1)

out1.close()
