#!/usr/bin/env python
#==============================================================================
# filter-matched-ref-mers.py
#
# Shawn Driscoll
# 20160211
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# After running make-junc-db.py, kmering the reference and matching the 
# kmers back to the reference with seal, this script is used to filter mers
# to drop those out that are uninformative and/or misleading. 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# load GTF annotation
	sys.stderr.write("loading juncdb GTF annotation...\n")
	annot = parse_gtf(args.gtf)
	sys.stderr.write("found {} features\n".format(len(annot.keys())))

	read_idx = 0
	mer_idx = 0
	mer_id = ""
	keep_mer = False
	tmp_type = ""

	fin = open(args.mers, "r")
	for szl in fin:
		szl = szl.strip()

		if szl[0]==">":
			read_idx += 1
			keep_mer = True
			aln = szl.split("\t")

			tnames = []
			thits = []
			types = set()
			gid = set()

			for i in range(1, len(aln)):
				tmp = aln[i].split("=")
				tnames.append(tmp[0])
				thits.append(int(tmp[1]))

				# deal with type. any mer that matches 5p and 3p and or a psi
				# site we don't want.
				tmp_type = annot[tmp[0]]['type']
				# if it is a theta then append 3p or 5p
				if tmp_type == "theta":
					res = re.search("([35]p)", tmp[0])
					if res:
						tmp_type += ":{}".format(res.group(1))

				types.update([tmp_type])
				gid.update([annot[tmp[0]]['gid']])

			if max(thits) > 1:
				keep_mer = False

			if len(types) > 1:
				keep_mer = False

			if len(gid) > 1:
				keep_mer = False

			if keep_mer:
				# revise this line so the name makes more sense
				mer_id = "MERID_{:08d}".format(mer_idx)
				mer_idx += 1
				target_set = "|".join(list(set(tnames)))
				szl = ">{}={}".format(mer_id, target_set)

		if keep_mer:
			print szl

	fin.close()
	

	sys.stderr.write("retained {} of {} mers\n".format(mer_idx, read_idx))

	return 0


#
# parser the GTF
def parse_gtf(fname):
	# variables
	gtfdb = {}
	szl = ""
	aln = []

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		# parse attributes
		attr = parse_gtf_attr(aln[8])

		if attr['transcript_id'] not in gtfdb:
			gtfdb[attr['transcript_id']] = {'gid': attr['gene_id'], 'gname': attr['gene_name'], 
			'tid': delim_list_to_set(attr['oId'], ","), 'hits': 0, 'type': aln[1], 'strand': aln[6]}

	fin.close()

	return gtfdb


def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	#
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

def delim_list_to_set(sz, delim):
	aa = sz.split(delim)
	da = {}

	for n in aa:
		da[n] = 0

	return da

# like the R order function which returns a sorting of the indices
def order(v):
	return np.argsort(v)


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('gtf', type=str, help="Reference jundb GTF")
parser.add_argument('mers', type=str, help="Reference kmers matched back against the reference [fasta format]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

