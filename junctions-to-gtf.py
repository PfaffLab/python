#!/usr/bin/env python
#==============================================================================
# junctions-to-gtf.py
#
# Shawn Driscoll
# 20120712
#
# Parses a junctions.bed file and produces a gtf where each junction
# is treated as a single isoform gene
#==============================================================================

import sys,argparse

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Parses a junctions.bed file into a GTF file treating each junction as a single isoform gene.")
parser.add_argument('bed',type=str,help="junctions BED file or enter - to read from stdin")

args = parser.parse_args()

#==============================================================================
# main script
#==============================================================================

arExon = [0,0]
arRow = ["","tophat_junctions","exon",0,0,0.,"",".",""]

if args.bed == "-":
	fin = sys.stdin
else:
	try:
		fin = open(args.bed,"r")
	except IOError as e:
		print "Error: cannot open input file: " + args.bed
		sys.exit()


# read from stdin
for sl in fin:
	sl = sl.strip()

	arl = sl.split("\t")

	if len(arl) > 1:
		#print len(arl)
		arExt = arl[10].split(",")
		arExt = map(int,arExt)
		
		arRow[0] = arl[0]
		arRow[3] = int(arl[1])
		arRow[4] = int(arl[1])+int(arExt[0])
		arRow[6] = "+"
		arRow[8] = "gene_id \"" + arl[3] + "\"; transcript_id \"" + arl[3] + "\"; FPKM \"" + arl[4] + "\";"

		print "\t".join(map(str,arRow))

		arRow[0] = arl[0]
		arRow[3] = int(arl[2])-arExt[1]
		arRow[4] = int(arl[2])

		print "\t".join(map(str,arRow))

if args.bed != "-":
	fin.close()

