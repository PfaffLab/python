#!/usr/bin/env python
#==============================================================================
# fasta-to-tab-delim.py
#
# Shawn Driscoll
# 20130201
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# 
#==============================================================================

import sys, argparse

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Generate a tab-delimited copy of your FASTA file.")
parser.add_argument('input_fa',type=str,help="FASTA format input file or - to read from STDIN")

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

	header = []

	#
	# read through file parsing out the reference names and their sequence
	fin = open(args.input_fa,"r")

	for szl in fin:
		szl = szl.strip()

		if szl[0] == ">":

			if len(header) > 0:
				print parse_entry(header, sz_seq)

			# reference line, if there is sequence loaded up then print out the 
			# GTF line
			header = szl[1:].split(" ")
			sz_seq = ""

		else:
			sz_seq += szl

	# done, close file
	fin.close()

	if len(header) > 0:
		print parse_entry(header, sz_seq)

	return 0

#==============================================================================
# makeGtfLine
# returns a string in GTF format based on supplied values
#==============================================================================
def parse_entry(header, seq):

	# let's see how much we can break the header up past just space delim
	new_line = []

	for i in range(len(header)):
		temp = header[i].split("=")
		new_line += temp

	new_line.append(seq)

	return "\t".join(new_line)



#==============================================================================
# entry point
#==============================================================================
if __name__ == "__main__":
	sys.exit(main(args))
