#!/usr/bin/env python
#==============================================================================
# fasta-to-gtf.py
#
# Shawn Driscoll
# 20120816
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Generate a GTF format annotation from a fasta file so that it may be used
# with a read counting script to quantify expression after alignment
#==============================================================================

import sys,argparse

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Generate a GTF annotation from a FASTA file - one entry per entry.")
parser.add_argument('input_fa',type=str,help="FASTA format input file or - to read from STDIN")
parser.add_argument("--bed",dest="b_bedout",action="store_const",const=True,default=False,help="Output BED format instead of GTF format (default: False)")

args = parser.parse_args()

# check for file
try:
	fin = open(args.input_fa,"r")
except IOError,e:
	sys.stderr.write("I/O Error ({0}): {1}\n".format(e.errno,e.strerror))
	sys.exit()

#==============================================================================
# globals
#==============================================================================

#==============================================================================
# main
#==============================================================================

def main(args):

	#
	# variables
	szl = ""
	sz_seq = ""
	fin = None
	sz_name = ""
	sz_chrom = ""
	sz_id = ""


	#
	# read through file parsing out the reference names and their sequence
	fin = open(args.input_fa,"r")

	for szl in fin:
		szl = szl.strip()

		if szl[0] == ">":
			# reference line, if there is sequence loaded up then print out the 
			# GTF line
			if len(sz_seq) > 0:
				if args.b_bedout:
					print makeBedLine(sz_chrom,sz_name,sz_id,len(sz_seq),"+")
				else:
					print makeGtfLine(sz_chrom,sz_name,sz_id,len(sz_seq),"+")

			sz_seq = ""
			sz_id,sz_name,sz_chrom = parseName(szl)

		else:
			sz_seq += szl

	# done, close file
	fin.close()

	# print final entry
	if len(sz_seq) > 0:
		if args.b_bedout:
			print makeBedLine(sz_chrom,sz_name,sz_id,len(sz_seq),"+")
		else:
			print makeGtfLine(sz_chrom,sz_name,sz_id,len(sz_seq),"+")

	return 0

#==============================================================================
# makeGtfLine
# returns a string in GTF format based on supplied values
#==============================================================================
def makeGtfLine(sz_chrom,sz_name,sz_id,n_length,sz_strand):
	#
	# variables
	l_out = []

	# fill l_out with output info
	l_out = [
		sz_chrom,
		"from_fasta",
		"gene",
		"1",
		"{:d}".format(n_length),
		"{:.6f}".format(0),
		sz_strand,
		".",
		"gene_id \"{}\"; transcript_id \"{}\";".format(sz_name,sz_id)]

	return "\t".join(l_out)


#==============================================================================
# makeGtfLine
# returns a string in GTF format based on supplied values
#==============================================================================
def makeBedLine(sz_chrom,sz_name,sz_id,n_length,sz_strand):
	#
	# variables
	l_out = []

	sz_nameCol = sz_id
	if sz_id != sz_name:
		sz_nameCol += "|" + sz_name

	# fill l_out with output info
	l_out = [
		sz_chrom,
		"0",
		"{:d}".format(n_length),
		sz_nameCol,
		"0",
		sz_strand]

	return "\t".join(l_out)

#==============================================================================
# parseName
# One might want to customize this function to properly parse different 
# fasta header names. By default text up to the first whitespace break is 
# returned for all three values. sz_chrom should always be set to whatever
# the reference name used by the aligner was. Best to check your alignment 
# files to see what was used.
#==============================================================================
def parseName(sz):
	#
	# variables
	sz_id = ""
	sz_name = ""
	sz_chrom = ""

	sz_id = sz[1:].split(" ")[0]
	sz_name = sz_id
	sz_chrom = sz_id
	return (sz_id,sz_name,sz_chrom)


#==============================================================================
# entry point
#==============================================================================
if __name__ == "__main__":
	sys.exit(main(args))
