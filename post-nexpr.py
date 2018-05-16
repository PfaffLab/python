#!/usr/bin/env python
#==============================================================================
# post-nexpr.py
#
# Shawn Driscoll
# 20120829
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Simplifies the process of merging several 'nexpr' generated expression files
# and generating the 'genes.counts.txt' and 'genes.expr.txt' files.
#==============================================================================

import sys,argparse

#==============================================================================
# globals
#==============================================================================

HITS_COL = -2
EXPR_COL = -1

#==============================================================================
# main
#==============================================================================

def main(args):
	
	#
	# variables
	counts_table = {}
	expr_table = {}
	file_index = 0
	num_files = 0
	header_names = []
	temp = []
	outfile_expr = "genes.expr.txt"
	outfile_counts = "genes.counts.txt"
	fp = None
	fp_a = None
	fp_b = None

		
	# passed, continue...
	num_files = len(args.expr_file)
	
	#
	# create header
	#
	
	header_names = ["tid","gid","location"]
	for i in range(num_files):
		temp = args.expr_file[i].split(".")
		header_names.append(temp[0])
	
	#
	# process files. if it's the first file then we need to build the tables
	#
	
	try:
		fp = open(args.expr_file[0],"r")
	except IOError,e:
		sys.stderr.write("failed to open input file ({:s})".format(args.expr_file[0]))
		return 1
	
	for szl in fp:
		szl = szl.strip()
		ll = szl.split("\t")
		
		if ll[args.link_col] not in counts_table:
			counts_table[ll[args.link_col]] = [ll[args.link_col],ll[args.link_col*-1+1],ll[2]]
			expr_table[ll[args.link_col]] = [ll[args.link_col],ll[args.link_col*-1+1],ll[2]]
		else:
			sys.stderr.write("warning: duplicate key found ({:s})\n".format(ll[args.link_col]))
		
		# if ll[LINK_COL] not in expr_table:
		#else:
		#	sys.stderr.write("warning: duplicate key found ({:s})\n".format(ll[LINK_COL]))
		
		counts_table[ll[args.link_col]].append(ll[HITS_COL])
		expr_table[ll[args.link_col]].append(ll[EXPR_COL])
	
	fp.close()
	
	for i in range(1,num_files):

		try:
			fp = open(args.expr_file[i],"r")
		except IOError,e:
			sys.stderr.write("failed to open input file ({:s})".format(args.expr_file[i]))
			return 1
		
		for szl in fp:
			szl = szl.strip()
			ll = szl.split("\t")
			
			if ll[args.link_col] not in counts_table:
				sys.stderr.write("warning: missing key found ({:s})".format(ll[args.link_col]))
			
			if ll[args.link_col] not in expr_table:
				sys.stderr.write("warning: missing key found ({:s})".format(ll[args.link_col]))
			
			counts_table[ll[args.link_col]].append(ll[HITS_COL])
			expr_table[ll[args.link_col]].append(ll[EXPR_COL])
		
		fp.close()
	
	#
	# finished merging data, now we can write the results to files
	#
	
	try:
		fp_a = open(outfile_counts,"w")
		fp_b = open(outfile_expr,"w")
	except IOError,e:
		sys.stderr.write("failed to open output files!\n")
		return 1
	
	#
	# write headers
	fp_a.write("\t".join(header_names) + "\n")
	fp_b.write("\t".join(header_names) + "\n")
	
	for key in sorted(counts_table.keys()):
		
		fp_a.write("\t".join(counts_table[key]) + "\n")
		fp_b.write("\t".join(expr_table[key]) + "\n")
		
	fp_a.close()
	fp_b.close()
	
	return 0

#==============================================================================
# additional defs
#==============================================================================

def check_args(args):
	
	if len(args) < 1:
		sys.stderr.write("too few arguments!\n")
		sys.exit(1)
		
	return 0

#==============================================================================
# entry point
#==============================================================================

#
# parse arguments
parser = argparse.ArgumentParser(description="Simplify merging of expression files output from nexpr")
parser.add_argument("expr_file",type=str,nargs="+",help="File(s) to parse")
parser.add_argument("-c",type=int,dest="link_col",default=1,help="Column index of unique keys, zero based (default: 1)")
args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))


