#!/usr/bin/env python
#==============================================================================
# remove-junction-aligns.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys, argparse, re

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	szl = ""
	ll = []
	rbuffer = []
	qname = ""
	qname_last = ""

	# check input file
	if args.infile != "-":
		if not file_exists(args.infile):
			sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.infile)
			return 1
	else:
		sys.stderr.write("[main] Expecting SAM input on STDIN...\n")


	# open dropped alignment file for reads where only junction alignments are 
	# reported
	fout = open("dropped.sam", "w")

	if args.infile == "-":
		fin = sys.stdin
	else:
		fin = open(args.infile, "r")

	for szl in fin:
		szl = szl.strip()

		if szl[0] == "@":
			print szl
		else:
			ll = szl.split("\t")
			qname = ll[0]

			if qname != qname_last:
				if len(rbuffer) > 0:
					# check buffer for juntion alignments

					kills = [0 for i in range(len(rbuffer))]

					for i in range(len(rbuffer)):
						cigar = rbuffer[i][5]

						if cigar.find("N") > 0:
							kills[i] = 1
					
					num_kills = sum(kills)
					if num_kills < len(rbuffer):
						# print out alignments that don't have to be killed
						for i in range(len(rbuffer)):
							if kills[i] == 0:
								sys.stdout.write("\t".join(rbuffer[i]) + "\n")

					else:
						# all alignments are bunk - send them all out to dropped.sam
						for i in range(len(rbuffer)):
							fout.write("\t".join(rbuffer[i]) + "\n")


				# clear buffer
				rbuffer = []

			rbuffer.append(ll)
			qname_last = qname

	# process final buffer

	if len(rbuffer) > 0:
		# check buffer for juntion alignments

		kills = [0 for i in range(len(rbuffer))]

		for i in range(len(rbuffer)):
			cigar = rbuffer[i][5]

			if cigar.find("N") > 0:
				kills[i] = 1
		
		num_kills = sum(kills)
		if num_kills < len(rbuffer):
			# print out alignments that don't have to be killed
			for i in range(len(rbuffer)):
				if kills[i] == 0:
					sys.stdout.write("\t".join(rbuffer[i]) + "\n")

		else:
			# all alignments are bunk - send them all out to dropped.sam
			for i in range(len(rbuffer)):
				fout.write("\t".join(rbuffer[i]) + "\n")


	fout.close()
	fin.close()

	return 0


def set_bit(fflag, value):
	new_flag = fflag | value
	return new_flag

def clear_bit(fflag, value):
	new_flag = fflag & ~value
	return new_flag

def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('infile', type=str, help="Input file")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
