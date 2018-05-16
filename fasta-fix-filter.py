#!/usr/bin/env python
#==============================================================================
# fasta-fix-filter.py
#
# Shawn Driscoll
# 20121127
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Reformats/filters a FASTA file. Ensures that all lines of the FASTA are the
# same width and allows the user to filter by total sequence length per
# feature.
#==============================================================================

import sys,argparse
import re

#==============================================================================
# main
#==============================================================================

def main(args):

	#
	# check file
	#
	try:
		fin = open(args.fasta_file,"r")
	except IOError,e:
		sys.stderr.write("Error: unable to open input file {:s}\n".format(args.fasta_file))
		return 1

	seq_name = ""
	seq = ""

	npad = ""
	if args.N > 0:
		tmp = ["N" for i in range(args.N)]
		npad = "".join(tmp)

	#
	# read and process the file
	#
	for szl in fin:
		# strip whitespace
		szl = szl.strip()
		if len(szl) < 1:
			continue
					
		# is this a sequence id line?
		if szl[0] == ">":
			if len(seq) > 0 and len(seq) > args.min_length:
				print_seq(seq_name, npad + seq + npad)

			seq_name = szl
			seq = ""
		else:
			if args.u:
				tmp = szl.upper()
				szl = tmp
			
			seq = seq + szl

	if len(seq) > 0 and len(seq) > args.min_length:
		print_seq(seq_name, npad + seq + npad)

	fin.close()

	return 0

def print_seq(name, seq):
	
	line_length = 50
	chunks = len(seq)/line_length

	if re.search("length\=[0-9]+", name):
		tmp = re.sub("[\s]length\=[0-9]+", "", name)
		name = tmp

	print name + " length=" + str(len(seq))

	for i in range(chunks):
		print seq[(i*line_length):((i+1)*line_length)]

	if len(seq) % line_length != 0:
		print seq[(line_length*chunks):]

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Reformats/filters a FASTA file. Ensures that all lines of the FASTA are the same length and provides the option to filter out sequences shorter than some cutoff value.")
parser.add_argument('fasta_file',type=str,help="FASTA format input file")
parser.add_argument('-f',dest="min_length",type=int,default=0,help="Optional sequence length cutoff. Only those sequences longer than this length are returned. (default: 0)")
parser.add_argument('-u', action="store_const", const=True, default=False, 
	help="Convert characters to upper case")
parser.add_argument('-N', type=int, default=0, 
	help="N-pad the sequences with this many N's at the start and end [0]")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
