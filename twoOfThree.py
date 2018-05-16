#!/usr/bin/env python
#==============================================================================
# twoOfThree
#
# shawn driscoll
# 20120710
#
# gene expression laboratory, pfaff
# salk institute for biological studies
#
# prints out rows from *.tracking files produced by cuffcompare
#==============================================================================

import sys,argparse,math

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Print rows from *.tracking file that meet some condition of number of samples with data (non - entries)")

parser.add_argument('infile',action='store',type=str,help="Parsed tracking file (enter - to read from stdin)")
parser.add_argument('--require-all',action='store_const',dest='bReqAll',const=True,default=False,help="Require all fields to contain a non-zero value (default:False)")
parser.add_argument('--require-none',action='store_const',dest='bReqNone',const=True,default=False,help="Require all fields to not have a value (default:False)")
parser.add_argument('-r',action='store',dest='nReq',default=-1,type=int,help="Require a specific number of samples to contain a value (default: majority)")
parser.add_argument('--parsed',action='store_const',dest='bParsed',default=False,const=True,help="Input tracking data has already been parsed by parse-tracking")

args = parser.parse_args()


#==============================================================================
# variables
#==============================================================================

szl = ""
arl = []
nScore = 0
nSamples = 0
nMajority = 0
nExpr = 0
nOffset = 4

#==============================================================================
# main script
#==============================================================================

if args.bParsed:
	nOffset = 5

# open input
if args.infile == "-":
	fin = sys.stdin
else:
	try:
		fin = open(args.infile,"r")
	except IOError as e:
		print "Error: cannot open input file: " + args.infile
		sys.exit()


for szl in fin:
	szl = szl.strip()
	arl = szl.split("\t")

	nSamples = len(arl)-nOffset
	if args.nReq == -1:
		nMajority = math.ceil(nSamples/2.)
	else:
		nMajority = args.nReq

	nExpr = 0

	nScore = nSamples
	for i in range(nSamples):
		if args.bParsed:
			if arl[i+nOffset] == "-" or float(arl[i+nOffset]) == 0:
				nScore -= 1
		else:
			if arl[i+nOffset] == "-":
				nScore -= 1


	if args.bReqAll:
		if nScore == nSamples:
			print szl
	elif args.bReqNone:
		if nScore == 0:
			print szl
	else:
		#print nMajority
		if nScore >= nMajority:
			print szl

if args.infile != "-":
	fin.close()



