#!/usr/bin/python
#==============================================================================
# fasta-index.py
#
# shawn driscoll
# 20120731
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Provide a fasta file and this script creates a special index file of 
# file pointer locations for each fasta entry in the file. Useful for indexing
# full genomes.
#==============================================================================

import sys,argparse,re

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Index a FASTA file for quick access.")
parser.add_argument('infile',type=str,help="FASTA index file")

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

n_plast = 0
l_row = ["",0,0,0]
n_lineChars = 0
n_lineActual = 0

sz_out = args.infile + "is"

#==============================================================================
# main script
#==============================================================================

fin = open(args.infile,"r")
fout = open(sz_out,"w")

# figure out line lengths
szl = fin.readline()
# second line has characters
szl = fin.readline()
l_row[2] = len(szl.strip())
l_row[3] = len(szl)
# reset to start of file
fin.seek(0)

# starting from top of file we want to index the file pointer positions of 
# the start and end of each reference's sequence

n_plast = 0
szl = fin.readline()
while szl:
	# look for fasta headers

	m = re.search("^>(.+)$",szl)
	if m:
		# we are at a header. if this is the first one then we just need to 
		# record the position of the file pointer because it is set to the start
		# of the current entry's sequence. if the last pointer isn't 0 then
		# we are at the end of a record and n_plast should contain the position
		# of the start of this row which is also the end of the last row plus
		# a character or two
		if n_plast > 0:
			# there was a previous reference - record its length
			#l_row[2] = fin.tell()-l_row[1]
			fout.write("\t".join(map(str,l_row)) + "\n")

		# record current reference name
		l_row[0] = m.group(1)
		# record start position of this reference's sequence
		l_row[1] = fin.tell()

	# keep current pointer before reading next line
	n_plast = fin.tell()
	szl = fin.readline()

if n_plast > 0:
	# there was a previous reference - record its length
	#l_row[2] = fin.tell()-l_row[1]
	fout.write("\t".join(map(str,l_row)) + "\n")

fin.close()
fout.close()

