#!/usr/bin/python
#==============================================================================
# kmer-mappability-to-graph.py
#
# Shawn Driscoll
# 20170316
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Matches the kmers of each fasta sequence in a fasta file to all 
# sequences in the file. Generates a count of the distinct kmers per 
# sequence as well as counts of shared kmers to other sequences.  
#==============================================================================

import sys, argparse
from math import log, exp
from os.path import isfile, expanduser
from os import system
# import re
# from igraph import *
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

	#
	# variables
	dpairs = {}
	nparsed = 0
	npairs = 0
	nkept = 0
	transcript_mode = False
	last_pairs = 0

	tmode_header = "tid1\ttid2\tsim\tshared_k\tinfo1\tinfo2\tintra_gene\tintra_gid"
	gmode_header = "gene1\tgene2\tsim\tshared_k\tinfo1\tinfo2\tintra_gene\tintra_gid"

	#
	# parse input file
	fin = open(args.mapp, "r")
	# skip header
	szl = fin.readline()
	aln = szl.strip().split("\t")

	# detect mode
	if aln[0] == "transcript_id":
		message("Transcript mappability detected")
		transcript_mode = True
	elif aln[0] == "group_id":
		message("Gene mappability detected")
		transcript_mode = False
	else:
		message("Cannot detect type of input file. Expected either 'transcript_id' or 'group_id' as the first column header.")
		return(1)

	if transcript_mode:
		print tmode_header
	else:
		print gmode_header
	
	for szl in fin:
		
		# parse
		if transcript_mode:
			trow = Tmapp(szl)

			if trow.num_shared > 0:
				# setup the pairs
				for htid in (trow.hits).keys():
					r = RefPair(trow.tid, htid)
					if r.key not in dpairs:
						# first encounter
						r.k1 = trow.klen
						r.shared = trow.hits[htid]
						r.gname1 = trow.gname
						r.gid1 = trow.gid
						dpairs[r.key] = r
					else:
						# second encounter
						dpairs[r.key].k2 = trow.klen
						dpairs[r.key].gname2 = trow.gname
						dpairs[r.key].gid2 = trow.gid
						dpairs[r.key].calc_sim()
						npairs += 1

						# pass count thresholds?
						print_pair = True
						if not dpairs[r.key].check(args.min_count, args.min_ratio):
							# nope
							print_pair = False

						if args.no_intra_gene & dpairs[r.key].is_intra_gene():
							print_pair = False
						
						if args.no_inter_gene & (not dpairs[r.key].is_intra_gene()):
							print_pair = False
						
						if args.no_intra_gid & dpairs[r.key].is_intra_gid():
							print_pair = False

						if args.no_inter_gid & (not dpairs[r.key].is_intra_gid()):
							print_pair = False
						
						if print_pair:
							nkept += 1
							print dpairs[r.key]

					if npairs != last_pairs:
						if (npairs % 1e5) == 0:
							last_paired = npairs
							message("parsed {} rows; found {} pairs; kept {} ({:0.1f}%)".format(nparsed, npairs, nkept, nkept*100.0/npairs))

		else:
			# gene mode
			grow = Gmapp(szl)

			if grow.num_shared > 0:
				# setup the pairs
				for hgid in (grow.hits).keys():
					r = RefPair(grow.group_id, hgid)
					if r.key not in dpairs:
						# first encounter
						r.k1 = grow.klen
						r.shared = grow.hits[hgid]
						r.gname1 = grow.gene_name
						r.gid1 = grow.gene_id
						dpairs[r.key] = r
					else:
						# second encounter
						dpairs[r.key].k2 = grow.klen
						dpairs[r.key].gname2 = grow.gene_name
						dpairs[r.key].gid2 = grow.gene_id
						dpairs[r.key].calc_sim()
						npairs += 1

						# pass count thresholds?
						print_pair = True
						if not dpairs[r.key].check(args.min_count, args.min_ratio):
							# nope
							print_pair = False

						if args.no_intra_gene & dpairs[r.key].is_intra_gene():
							print_pair = False
						
						if args.no_inter_gene & (not dpairs[r.key].is_intra_gene()):
							print_pair = False
						
						if args.no_intra_gid & dpairs[r.key].is_intra_gid():
							print_pair = False

						if args.no_inter_gid & (not dpairs[r.key].is_intra_gid()):
							print_pair = False
						
						if print_pair:
							nkept += 1
							print dpairs[r.key]

					if npairs != last_pairs:
						if (npairs % 1e5) == 0:
							last_pairs = npairs
							message("parsed {} rows; found {} pairs; kept {} ({:0.1f}%)".format(nparsed, npairs, nkept, nkept*100.0/npairs))


		nparsed += 1

	fin.close()

	message("parsed {} rows; found {} pairs; kept {} ({:0.1f}%)".format(nparsed, npairs, nkept, nkept*100.0/npairs))

	return 0


#==============================================================================
# defs
#==============================================================================

#
# class for transcript mappability row
class Tmapp(object):

	def __init__(self, sz):
		
		aln = sz.strip().split("\t")
		
		self.tid = aln[0]
		self.gid = aln[1]
		self.gname = aln[2]
		self.klen = int(aln[8])

		self.hits = {}

		if int(aln[9]) > 0:
			tmp = aln[11].split(";")
			tmp2 = aln[12].split(";")
			for i in range(len(tmp)):
				self.hits[tmp[i]] = int(tmp2[i])

		if int(aln[10]) > 0:
			tmp = aln[13].split(";")
			tmp2 = aln[14].split(";")
			for i in range(len(tmp)):
				self.hits[tmp[i]] = int(tmp2[i])

		self.num_shared = len((self.hits).keys())

		return None


#
# class for gene type mappability row
class Gmapp(object):

	def __init__(self, sz):
		
		aln = sz.strip().split("\t")
		
		self.group_id = aln[0]
		self.gene_id = aln[1]
		self.gene_name = aln[2]
		self.transcripts = aln[3]
		self.mapp = float(aln[4])
		self.min_mapp = float(aln[5])
		self.klen = int(aln[6])

		self.hits = {}

		if int(aln[7]) > 0:
			tmp = aln[8].split(";")
			tmp2 = aln[9].split(";")
			for i in range(len(tmp)):
				self.hits[tmp[i]] = int(tmp2[i])

		self.num_shared = len((self.hits).keys())

		return None


#
# class to hold pair of reference terms, their kmer lengths and their shared count
class RefPair(object):

	def __init__(self, id1, id2):
		# init class attributes
		self.id1 = id1
		self.id2 = id2

		self.k1 = 0
		self.k2 = 0
		self.shared = 0

		self.gname1 = "u"
		self.gname2 = "u"
		self.gid1 = "u"
		self.gid2 = "u"

		tmp = sorted([id1, id2])
		self.key = ";".join(tmp)

		self.sim = 0

		return None

	def __str__(self):
		ltmp = self.tolist()
		return "\t".join(map(str, ltmp))

	def tolist(self):
		r1, r2 = self.calc_ratios()
		info1 = "{}|{}|{}|{:0.4f}".format(self.gid1, self.gname1, self.k1, r1)
		info2 = "{}|{}|{}|{:0.4f}".format(self.gid2, self.gname2, self.k2, r2)

		gname_match = 0
		if self.is_intra_gene():
			gname_match = 1
		gid_match = 0
		if self.is_intra_gid():
			gid_match = 1

		lout = [self.id1, self.id2, self.sim, self.shared, 
			info1, info2, gname_match, gid_match]
		return(lout)

	def is_intra_gene(self):
		return self.gname1==self.gname2

	def is_intra_gid(self):
		return self.gid1==self.gid2

	def calc_ratios(self):
		r1 = self.shared*1.0/self.k1
		r2 = self.shared*1.0/self.k2
		return r1, r2		

	def calc_sim(self):
		r1, r2 = self.calc_ratios()
		self.sim = pair_sim(r1, r2, self.k1, self.k2)
		return 0

	# compare the values in this pair to a minimum shared kmer 
	# level as well as a minimum shared ratio
	def check(self, min_n, min_r):
		rres = False
		r1, r2 = self.calc_ratios()
		
		if (self.shared >= min_n) and (r1 > min_r or r2 > min_r):
			rres = True

		return rres


#
# weighted mean
def wmean(v, w):

	n = len(v)
	if n != len(w):
		return 0

	nsum = 0
	ndenom = 0
	for i in range(n):
		nsum += v[i]*w[i]
		ndenom += w[i]

	return nsum/ndenom



def mean(v):

	n = len(v)
	vsum = 0
	vmean = 0

	for i in range(n):
		vsum += v[i]

	vmean = vsum/float(n)

	return(vmean)


#
# given the sim scores between two transcripts we can summarize it into a 
# single score with this function by producing a weighted mean based on the 
# lengths of the transcripts. 
# parameters:
# s1 is ratio of t1 shared by t2
# s2 is ratio of t2 shared by t1
def pair_sim(s1, s2, l1, l2):
	#
	# create weights from the lengths
	ll1 = log(l1)
	ll2 = log(l2)
	llmean = mean([ll1, ll2])
	w1 = 1/exp(ll1-llmean)
	w2 = 1/exp(ll2-llmean)

	rres = wmean([s1, s2], [w1, w2])

	return rres


def ftos(x):
	return("{:0.4f}".format(x))

def quantile(v, q):
	n = len(v)
	vhat = sorted(v)

#
# run a system level command. used for running alignment and samtools commands
def runcmd(cmd, verbose=True):
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	rres = system(cmd)
	return rres

#
# print a message to STDERR. appends the newline for you.
def message(sz):
	sys.stderr.write("[kmer-gene-mappability] " + sz + "\n")


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.")

# -- required args
parser.add_argument('mapp', type=str, help="Mappability output")

# -- options
parser.add_argument('-c', '--min-count', type=int, default=1, 
	help="Minimum number of kmers to keep as a match. Used in combination with -r. [1]")
parser.add_argument('-r', '--min-ratio', type=float, default=0.001, 
	help="Minimum ratio of kmers shared to keep as a match. Used in combination with -c. [0.001]")
parser.add_argument('--no-intra-gene', action="store_const", const=True, default=False, 
	help="Filter output. Do not print intra-gene pairs.")
parser.add_argument('--no-intra-gid', action="store_const", const=True, default=False, 
	help="Filter output. Do not print intra-gene_id pairs.")
parser.add_argument('--no-inter-gene', action="store_const", const=True, default=False, 
	help="Filter output. Do not print inter-gene pairs.")
parser.add_argument('--no-inter-gid', action="store_const", const=True, default=False, 
	help="Filter output. Do not print inter-gene_id pairs.")

# -- parse command line
args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

