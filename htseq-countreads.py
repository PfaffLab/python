#!/usr/bin/env python
#==============================================================================
# fastq-countReads.py
#
# Shawn Driscoll
# 20120502
#
# counts reads collected in FASTQ file (may come from stdin)
#==============================================================================

import sys,argparse
import HTSeq
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# parse args
#==============================================================================

parser = argparse.ArgumentParser(description="Count reads in a FASTQ file")
parser.add_argument("infile",action="store",type=str,help="FASTQ input or - to read from stdin")
args = parser.parse_args()

#==============================================================================
# main scirpt
#==============================================================================

if args.infile == "-":
	in1 = iter(HTSeq.FastqReader(sys.stdin))
else:
	in1 = iter(HTSeq.FastqReader(args.infile))

nReads = 0
nReadLength = 0

for read in in1:
	nReads += 1

	if nReadLength == 0:
		nReadLength = len(read)
		arQuals = np.zeros(len(read),np.int)

	arQuals += read.qual


# compute mean qualities
arQuals /= float(nReads)

# plot qualities as a line plot
plt.plot(arQuals)
plt.ylabel("Base quality score")
plt.xlabel("Base position in reads")
plt.savefig("fastq-qualitiesPlot.pdf",dpi=300,format="pdf")

print "> reads:       " + str(nReads)
print "> read length: " + str(nReadLength) + "bp"
