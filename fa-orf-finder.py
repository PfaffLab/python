#!/usr/bin/env python
#==============================================================================
# fa-orf-finder.py
#
# shawn driscoll
# 20120711
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Expects a FASTA format file. Searches for longest open reading frame within
# each record.
#==============================================================================

import sys, argparse, re

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(
	description="Searches for longest open reading frame within each FASTA record. Returns FASTA format file with either the ORF sequenc found or a short sequence of N's")
parser.add_argument('infile',type=str,help="FASTA format file")
parser.add_argument('--fa-cds-only',dest='b_cdsOnly',action="store_const",const=True,default=True,
	help="Return CDS sequence only otherwise return original sequence (default: True)")

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

bedout = args.infile + ".orf.bed"
faout = args.infile + ".orf.fa"
szSeq = ""
szHeader = ""
nStart = 0
nEnd = 0
nseq = 0

#==============================================================================
# functions
#==============================================================================

# -- 
# searches sequence for open reading frames, returns starts and ends
# of each one found (single sequence may have many)
def orfAnalysis(szSeq, szHeader, b_cdsOnly, f):
	nLen = len(szSeq)
	szSeq = szSeq.upper()
	szStart = "ATG"
	dStops = {"TAG":0,"TGA":0,"TAA":0}
	arFoundStops = [-1,-1,-1]
	bNoStart = True
	bOrf = False
	bFound = True
	ni = 0
	nj = 0

	arStarts = []
	arEnds = []
	arLengths = []

	# find all start codons
	ni = 0
	bFound = True
	while ni < nLen-2:
		if szSeq[ni:(ni+3)] == szStart:
			arStarts.append(ni)
			ni += 3
		else:
			ni += 1

	if len(arStarts) == 0:
		# no starts so we can drop out right here
#		f.write(szHeader + " CDS=-1\n")
#		f.write(szSeq + "\n")
		return ([-1],[-1])

	# initalize lists for the ends and lenghts of any ORF starting from
	# the already discovered starts
	arEnds = [-1 for x in range(len(arStarts))]
	arLengths = [-1 for x in range(len(arStarts))]

	# for each start search for stop codons within frame
	for nj in range(len(arStarts)):
		# start at the next codon
		ni = arStarts[nj]+3
		bFound = False
		# loop through the sequence until we either find an end or we run out of 
		# sequence
		while ni < nLen-2 and not bFound:
			# check for all three stops
			if szSeq[ni:(ni+3)] in dStops:
				arEnds[nj] = ni
				bFound = True
			# go to next codon
			ni += 3


		if not bFound:
			# no end found for this start
			arEnds[nj] = -1
			arLengths[nj] = -1
		else:
			# found a stop codon for this start, find the length
			arLengths[nj] = arEnds[nj]-arStarts[nj]
			if arLengths[nj] % 3 != 0:
				# something is wrong...
				arEnds[nj] = -1
				arLengths[nj] = -1
				sys.stderr.write("Dropped ORF because length was not a multiple of 3\n")

	#print arStarts
	#print arEnds
	#print arLengths

	# figure out which ORF is longest
	nMax = arLengths[0]
	niMax = 0
	for ni in range(1,len(arStarts)):
		if arLengths[ni] > nMax:
			niMax = ni
			nMax = arLengths[ni]

	if nMax > 0:
		f.write(szHeader + " CDS=" + str(arStarts[niMax]+1) + "-" + str(arEnds[niMax]) + "\n")
		if b_cdsOnly:
			f.write(szSeq[arStarts[niMax]:arEnds[niMax]+3] + "\n")
		else:
			f.write(szSeq + "\n")
#	else:
#		f.write(szHeader + " CDS=-1\n")
#		f.write(szSeq + "\n")

	# RETURN lists of starts and ends
	return (arStarts,arEnds)


#==============================================================================
# main script
#==============================================================================

fin = open(args.infile,"r")

fout_fa = open(faout, "w")

fout_bed = open(bedout, "w")

for szl in fin:
	szl = szl.strip()

	# is this a seq header?
	if szl[0] == ">":
		nseq += 1

		if nseq % 1000 == 0:
			sys.stderr.write("Processed {:d} sequences\n".format(nseq))

		# start new seq
		if len(szSeq) > 0:
			# do orf analysis on previous seq
			lStart,lEnd = orfAnalysis(szSeq, ">" + szHeader, args.b_cdsOnly, fout_fa)
			lLens = [0 for i in range(len(lStart))]
			maxLen = 0
			maxI = -1
			for i in range(len(lStart)):
				if lEnd[i] > 0:
					lLens[i] = lEnd[i]-lStart[i]
					if lLens[i] > maxLen:
						maxLen = lLens[i]
						maxI = i

				lStart[i] += 1
			nStart = ",".join(map(str, lStart))
			nEnd = ",".join(map(str, lEnd))
			nLens = ",".join(map(str, lLens))
			fout_bed.write("\t".join(map(str,[szHeader, maxLen, maxI, nStart, nEnd, nLens])) + "\n")
		szSeq = ""
		htmp = re.sub("^>", "", szl)
		szHeader = re.search("(^[^\s]+)", htmp).group(1)

	else:
		szSeq += szl

fin.close()

# do orf analysis on final seq
lStart,lEnd = orfAnalysis(szSeq, ">" + szHeader, args.b_cdsOnly, fout_fa)
lLens = [0 for i in range(len(lStart))]
maxLen = 0
maxI = -1
for i in range(len(lStart)):
	if lEnd[i] > 0:
		lLens[i] = lEnd[i]-lStart[i]
		if lLens[i] > maxLen:
			maxLen = lLens[i]
			maxI = i

	lStart[i] += 1

nStart = ",".join(map(str, lStart))
nEnd = ",".join(map(str, lEnd))
nLens = ",".join(map(str, lLens))
fout_bed.write("\t".join(map(str,[szHeader, maxLen, maxI, nStart, nEnd, nLens])) + "\n")

fout_fa.close()
fout_bed.close()

sys.stderr.write("Processed {:d} sequences\n".format(nseq))
sys.stderr.write("Done!\n")
