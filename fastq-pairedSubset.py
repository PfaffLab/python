#!/usr/bin/env python
#==============================================================================
# fastq-pairedSubset.py
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

parser = argparse.ArgumentParser(description="Select randomly select reads from paired FASTQ files by percentage of total reads in file")
parser.add_argument("fraction",type=float,action="store",help="Fraction of reads as a decimal (range: 0,1)")
parser.add_argument("infile_1",action="store",type=str,help="FASTQ input 1")
parser.add_argument("infile_2",action="store",type=str,help="FASTQ input 2")
parser.add_argument("outfile_1",action="store",type=str,help="FASTQ output 1")
parser.add_argument("outfile_2",action="store",type=str,help="FASTQ output 2")
args = parser.parse_args()

#==============================================================================
# main scirpt
#==============================================================================

in1 = iter(HTSeq.FastqReader(args.infile_1))
in2 = iter(HTSeq.FastqReader(args.infile_2))
out1 = open(args.outfile_1, "w")
out2 = open(args.outfile_2, "w")

for read1,read2 in itertools.izip(in1,in2):
	if random.random() < args.fraction:
		read1.write_to_fastq_file(out1)
		read2.write_to_fastq_file(out2)

out1.close()
out2.close()
