#!/usr/bin/env python
#==============================================================================
# post-nexpr-exons.py
#
# Shawn Driscoll
# 20120829
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Simplifies the process of merging several 'nexpr-exons' generated expression files
# and generating the exons.counts.txt file.
#==============================================================================

import sys,argparse

#==============================================================================
# globals
#==============================================================================


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
	outfile_counts = "exons.counts.txt"
	fp = None
	fp_a = None
	fp_b = None
	tid = ""

	i = 0
	j = 0
	exon_hits = []
	exon_lengths = []

		
	# passed, continue...
	num_files = len(args.expr_file)
	
	#
	# create header
	#
	
	header_names = ["texon_id","gene_name","location"]
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

		tid = ll[args.link_col]
		exon_hits = ll[5].split(",")
		exon_lengths = ll[4].split(",")

		for i in range(len(exon_hits)):
			key = tid + ":{:03d}".format(i)

			if key in counts_table:
				sys.stderr.write("what the? duplicate key - how is that possible? ({:s})".format(key))

			counts_table[key] = [key,ll[1],ll[2],exon_hits[i]]

	
	fp.close()
	
	for j in range(1,num_files):

		try:
			fp = open(args.expr_file[j],"r")
		except IOError,e:
			sys.stderr.write("failed to open input file ({:s})".format(args.expr_file[j]))
			return 1
		
		for szl in fp:
			szl = szl.strip()
			ll = szl.split("\t")

			tid = ll[args.link_col]
			exon_hits = ll[5].split(",")
			exon_lengths = ll[4].split(",")

			for i in range(len(exon_hits)):
				key = tid + ":{:03d}".format(i)

				if key not in counts_table:
					sys.stderr.write("Error: missing key from file {:s} ({:s})".format(args.expr_file[i],key))

				counts_table[key].append(exon_hits[i])

		fp.close()
	
	#
	# finished merging data, now we can write the results to files
	#
	
	try:
		fp_a = open(outfile_counts,"w")
	except IOError,e:
		sys.stderr.write("failed to open output file!\n")
		return 1
	
	#
	# write headers
	fp_a.write("\t".join(header_names) + "\n")
	
	for key in sorted(counts_table.keys()):
		
		fp_a.write("\t".join(counts_table[key]) + "\n")
		
	fp_a.close()
	
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
parser.add_argument("-c",type=int,dest="link_col",default=0,help="Column index of transcript keys, zero based (default: 0)")
parser.add_argument("-f",type=int,dest="low_filter",default=0,help="Filter out rows where there are not at least 1/2 of the rows with this many counts. (default: 0)")
args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))


