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
import numpy as np

#==============================================================================
# parse args
#==============================================================================

parser = argparse.ArgumentParser(description="Count reads in a FASTQ file")
parser.add_argument("infile",action="store",type=str,help="FASTQ input or - to read from stdin")
parser.add_argument("--phred33",dest="phred33",action="store_const",const=True,default=True,help="qualities are Phred+33 (default)")
parser.add_argument("--phred64",dest="phred64",action="store_const",const=True,default=False,help="qualities are Phred+64")
parser.add_argument("--quals",dest="doQuals",action="store_const",const=True,default=False,help="produce base quality statistics (slow) (default: False)")
args = parser.parse_args()

#==============================================================================
# main scirpt
#==============================================================================

if args.infile == "-":
	fin = sys.stdin
else:
	fin = open(args.infile,"r")

nReads = 0.0
nReadLength = 0
nReadl = 0
ni = 0
nMeanRead = 0
arQuals = []
nQuals = 0.
lQualQuartiles = []
lQualValues = [0,0,0,0,0]

nQualScale = 33
if args.phred64:
	nQualScale = 64


for szl in fin:
	szl = szl.strip()
	ni += 1

	if ni==1:
		nReads += 1
	elif ni==4:
		if args.doQuals:
			arQuals = np.array([ord(szl[i])-nQualScale for i in range(len(szl))])
			nQuals += arQuals.mean()
			for i in range(5):
				lQualValues[i] += arQuals[lQualQuartiles[i]]

		ni=0

	elif ni==2:
		if args.doQuals:
			if nReadl==0:
				nReadl=len(szl)
				lQualQuartiles = [0,int(nReadl*0.25),int(nReadl*0.5)-1,int(nReadl*0.75)-1,nReadl-1]
		nReadLength += len(szl)


if args.infile != "-":
	fin.close()

# find mean read length
nMeanRead = nReadLength/nReads
# find mean base quality
if args.doQuals:
	for i in range(5):
		lQualValues[i] /= nReads

# print results

print "total reads (millions)\t" + str(nReads/1000000.)
print "total Mb read\t" + str(nReadLength/1000000.)
print "mean read length\t" + str(nMeanRead)
if args.doQuals:
	print "qual quartiles\t" + ",".join(map(str,map(int,lQualValues)))
