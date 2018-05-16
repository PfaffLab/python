#!/usr/bin/env python
#==============================================================================
# proto-ja-explode.py
#
# shawn driscoll
# 20120803
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# 
#==============================================================================

import sys,argparse,re
import HTSeq as hts
import pybedtools as pbr
import numpy as np
import hashlib

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Explode proto 'ja' file so it is organized by transcript id instead of junction id")
parser.add_argument('infile',type=str,help="proto.ja.txt input file")

args = parser.parse_args()

# check for file
try:
	fin = open(args.infile,"r")
	fin.close()
except IOError as e:
	print "Error: unable to find/open input file: " + args.infile
	sys.exit()

#==============================================================================
# variables
#==============================================================================

dTids = {}
szl = ""
ll = []
lTids = []
szCompact = ""

#==============================================================================
# main script
#==============================================================================

##
# read through the input file and build up the dTids dict with all of the transcript
# id's found in the file. at each transcript id we'll store a list of both the annotated junctions
# as well as the novel type junctions found for the transcript

fin = open(args.infile,"r")

for szl in fin:

	szl = szl.strip()
	ll = szl.split("\t")
	lTids = ll[3].split(",")

	for szTid in lTids:
		if szTid not in dTids:
			dTids[szTid] = {"annot":[], "novel":[]}

		# check type
		if ll[4] == "=":
			szCompact = "|".join([ll[0],ll[1],ll[5]])
			dTids[szTid]["annot"].append(szCompact)
		else:
			szCompact = "|".join([ll[0],ll[1],ll[5],ll[6],ll[7],ll[8]])
			dTids[szTid]["novel"].append(szCompact)


fin.close()

##
# print back out sorted by transcript id
for szTid in sorted(dTids.keys()):

	sys.stdout.write(szTid + "\t")
	if len(dTids[szTid]['annot']) > 0:
		sys.stdout.write(",".join(dTids[szTid]['annot']) + "\t")
	else:
		sys.stdout.write("-\t")

	if len(dTids[szTid]['novel']) > 0:
		sys.stdout.write(",".join(dTids[szTid]['novel']) + "\t")
	else:
		sys.stdout.write("-\t")
	
	sys.stdout.write("\n")
