#!/usr/bin/env python
#==============================================================================
# base-stats.py
#
# Shawn Driscoll
# 20120726
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Process the output of samtools 'mpileup' command and count base calls per 
# A/C/T/G at each base in the report.
#==============================================================================

import sys,argparse

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Process the output of samtools 'mpileup' command and count base calls per A/C/T/G at each base in the report.")
parser.add_argument('infile',type=str,help="Input file or - to read from stdin (pipe straight from mpileup)")

args = parser.parse_args()

# check for file
if args.infile != "-":
	try:
		fin = open(args.infile,"r")
	except IOError as e:
		print "Error: unable to find/open input file: " + args.infile
		sys.exit()

#==============================================================================
# variables
#==============================================================================

dBases = {}
sAccepted = set(['A','T','C','G','.',',','a','t','c','g'])
szb = ""
szbout = ""

#==============================================================================
# main script
#==============================================================================


if args.infile == "-":
	fin = sys.stdin
else:
	fin = open(args.infile,"r")


for szl in fin:
	szl = szl.strip()
	ll = szl.split("\t")

	# print first 4 columns back out
	sys.stdout.write("\t".join(ll[0:4]) + "\t")

	# process 5th column
	szb = ll[4]

	# loop through the string
	dBases = {}
	for i in range(len(szb)):

		if szb[i] in sAccepted:
			if szb[i] not in dBases:
				dBases[szb[i]] = 0

			dBases[szb[i]] += 1

	print dBases

	#szbout = ""
	#for key in sorted(dBases.keys()):
	#	if len(szbout) > 0:
	#		szbout += "\t"

	#	szbout += "\t".join(map(str,[key,dBases[key]]))

	#print szbout



if args.infile != "-":
	fin.close()

