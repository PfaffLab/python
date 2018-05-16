#!/usr/bin/env python
#==============================================================================
# mirna-filter-matched-ref-mers.py
#
# Shawn Driscoll
# 20160519
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# adapted from juncdb-filter-matched-ref-mers.py. 
# for building a micro-rna kmer index for mapping reads.
# 
# $ kmercountexact.sh in=index.fa out=mers.fa trd=t rcomp=f k=N
# $ seal.sh in=mers.fa ref=index.fa outm=mers_to_index_N.fa k=N hdist=0 mm=f rcomp=f rename=t
# $ mirna-filter-matched-ref-mers.py mers_to_index_N.fa > mer_index.N.fa
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

			for i in range(1, len(aln)):
				tmp = aln[i].split("=")
				tnames.append(tmp[0])
				thits.append(int(tmp[1]))

#			if max(thits) > 1:
#				keep_mer = False
#

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
parser.add_argument('mers', type=str, help="Reference kmers matched back against the reference [fasta format]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

