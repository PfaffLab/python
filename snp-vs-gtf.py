#!/usr/bin/env python
#==============================================================================
# snp-vs-gtf.py
#
# Shawn Driscoll
# 20130618
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Intersects a SNP file (like a VCF or output from exactSnp) with a GTF file.
#==============================================================================

import sys, argparse

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

	annotate = args.attr_list is not None
	lattrs = []
	attr_sets = []
	num_attrs = 0
	lheader = []
	szl = ""

	#----------- check input files

	if not file_exists(args.snp_file):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.snp_file)
		return 1
	if not file_exists(args.gtf_file):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.gtf_file)
		return 1

	#----------- get to work

	sys.stderr.write("[main] parsing GTF file...\n")
	lk_table = gtf_lookup_hash(args.gtf_file, args.feature_type)

	# open snp file and find stuff

	sys.stderr.write("[main] intersecting SNPs with GTF features...\n")

	fin = open(args.snp_file, "r")

	# deal with the header
	szl = fin.readline()
	lheader = szl.strip().split("\t")

	if annotate:
		lattrs = args.attr_list.split(",")
		num_attrs = len(lattrs)
		for attr in lattrs:
			lheader.append(attr)

	# print header 
	print "\t".join(lheader)

	# loop through snp rows
	for szl in fin:
		ll = szl.strip().split("\t")

		if len(ll) > 1:

			# find intersections
			lhits = lookup_position(lk_table, ll[0:2])
			
			if args.invert and len(lhits) == 0:
				# print this one back out
				print "\t".join(ll)
			elif not args.invert and len(lhits) > 0:
				# print this one out
				if annotate:
					# make sets for each attribute we want to keep
					attr_sets = [set([]) for i in range(num_attrs)]
					# populate the sets
					for i in range(len(lhits)):
						attr = parse_gtf_attr(lhits[i][8])

						for j in range(num_attrs):
							if lattrs[j] in attr:
								attr_sets[j].update([attr[lattrs[j]]])
							else:
								sys.stderr.write("[main] warning: gtf attr does not exist\n")

					for i in range(num_attrs):
						ll.append(",".join(list(attr_sets[i])))

				print "\t".join(ll)

	fin.close()

	return 0


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

# 
# parse_gtf_attr
def parse_gtf_attr(field):

	attrs = {}

	# split into fields
	l1 = field.split(";")

	# split each field into key and value
	l2 = []
	for i in range(len(l1)):
		l2.append(l1[i].split("\""))


	for i in range(len(l2)):
		if len(l2[i]) > 1:
			attrs[l2[i][0].strip()] = l2[i][1].strip()

	return(attrs)		

#
# make a fast lookup hash of the GTF
#
def gtf_lookup_hash(fname, ftype):

	lk = {}

	#- open file -#

	fin = open(fname, "r")

	#- parse file and build hash

	for szl in fin:
		ll = szl.strip().split("\t")
		
		if ll[2] == ftype:
			kid = hash_pos(ll[0], ll[3])

			if kid not in lk:
				lk[kid] = []

			lk[kid].append(list(ll))

	fin.close()

	return(lk)

def hash_pos(rname, pos):
	mod = 16000
	bin = int(pos)/mod
	kid = rname + str(bin)
	return kid

# 
# look up a position in the hash, return list of intersections
#
def lookup_position(lk, lpos):

	result = []

	query_id = hash_pos(lpos[0], lpos[1])

	if query_id in lk:
		# check for intersections with features in this bucket
		ll = lk[query_id]
		for i in range(len(ll)):
			lf = ll[i]
			if int(lf[3]) <= int(lpos[1]) and int(lf[4]) >= int(lpos[1]):
				# append entire GTF record so that the main function can 
				# have access to all feature information
				result.append(list(lf))

	return result


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Find hits (or non-hits) for SNPs relative to a GTF file")
parser.add_argument('snp_file', type=str, help="SNP input file")
parser.add_argument('gtf_file', type=str, help="GTF reference file")
parser.add_argument('-v', dest="invert", action="store_const", const=True, default=False, help="Return features that do not intersect anything in the GTF")
parser.add_argument('-a', type=str, dest='attr_list', action="store", default=None, help="Annotate the features with intersections with GTF attributes")
parser.add_argument('-f', type=str, dest='feature_type', action="store", default="exon", 
	help="Feature type to find intersections (default: exon)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
