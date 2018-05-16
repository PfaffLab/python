#!/usr/bin/python
#==============================================================================
# kmer-transcript-mappability.py
#
# Shawn Driscoll
# 20170315
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Matches the kmers of each fasta sequence in a fasta file to all 
# sequences in the file. Generates a count of the distinct kmers per 
# sequence as well as counts of shared kmers to other sequences.  
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
from os import system

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
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

	# variables
	szfai = "{}.fai".format(args.fasta)
	rres = None
	tid2gid = {}
	tid2gname = {}
	has_info = False
	fout_name = "{}.tmapp".format(args.fasta)
	ntid = 0

	#
	# check for fasta index
	if not fai_exists(args.fasta):
		# need to build on
		message("building FAI index")
		build_fai(args.fasta)

	#
	# parse the fai
	dfai = parse_fai(szfai)
	
	#
	# if info is supplied load it
	if args.info is not None:
		message("Loading annotation from {}".format(args.info))
		tid2gid, tid2gname = parse_info(args.info)
		has_info = True

	# 
	# this will hold the per-transcript results
	dtid = {}

	fout = open(fout_name, "w")

	#
	# test this idea with the first few
	all_tid = sorted(dfai.keys())
	n = len(all_tid)
	
	for tid in dfai.keys():

		if dfai[tid] >= args.kmer_length:

			#
			# fetch sequence and kmer it
			message("Processing {}".format(tid))
			cmd = "samtools faidx {} {} > tmp.fa".format(args.fasta, tid)
			runcmd(cmd, args.v)
			cmd = "kmercountexact.sh in=tmp.fa k={} rcomp=f out=tmp.mer.fa overwrite=t 2>kmer.log".format(args.kmer_length)
			runcmd(cmd, args.v)
			# now we have to match the kmers against the full fasta using seal
			cmd = "seal.sh in=tmp.mer.fa ref={} k={} threads={} rcomp=f mm=f hdist=0 out=tmp.hits.fa rename=t ambig=all trd=t overwrite=t 2>seal.log".format(args.fasta, args.kmer_length, args.threads)
			runcmd(cmd, args.v)
			
			message("parsing hits...")
			if has_info:
				rres = parse_ii_hits("tmp.hits.fa", tid, dtid, tid2gid, tid2gname, dfai)
			else:
				rres = parse_hits("tmp.hits.fa", tid, dtid, dfai)

			runcmd("rm tmp.fa tmp.mer.fa tmp.hits.fa", args.v)
		else:
			# can't process features that are shorter than the kmer length
			if has_info:
				rres = Tobj(tid, gid=tid2gid[tid], gname=tid2gname[tid], length=dfai[tid])
			else:
				rres = Tobj(tid, length=dfai[tid])

		ntid += 1
		if ntid==1:
			fout.write(rres.pheader() + "\n")

		lout = rres.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")

	fout.close()

	# print results
	#rres = print_results(dtid, fout_name)

	return 0


#==============================================================================
# defs
#==============================================================================

class Tobj(object):
	def __init__(self, tid, gid="u", gname="u", length=0):

		self.id = tid
		self.gid = gid
		self.gname = gname
		self.length = length
		self.num_k = 0

		self.inter_hits = {}
		self.intra_hits = {}
		
		self.num_inter = 0
		self.num_intra = 0


	def __str__(self):
		lout = self.tolist()
		return "\t".join(map(str, lout))

	def add_inter_hit(self, tid):
		if tid not in self.inter_hits:
			self.inter_hits[tid] = 0
			self.num_inter += 1

		self.inter_hits[tid] += 1
		return 0

	def add_intra_hit(self, tid):
		if tid not in self.intra_hits:
			self.intra_hits[tid] = 0
			self.num_intra += 1

		self.intra_hits[tid] += 1
		return 0

	def same_gid(self, gid):
		return(gid==self.gid)

	# increment the kmer counter
	def kpp(self):
		self.num_k += 1

	def hits_to_lists(self, type="inter"):

		foo = True
		hlist = []
		hvals = []
		hrats = []

		if type=="inter":
			if self.num_inter == 0:
				foo = False
			hits = self.inter_hits
		else:
			if self.num_intra == 0:
				foo = False
			hits = self.intra_hits

		if foo:
			hlist = sorted(hits.keys())
			hvals = []
			hrats = []
			for k in hlist:
				hvals.append(hits[k])
				hrats.append(hits[k]/float(self.num_k))

		return hlist, hvals, hrats

	# calculate mappability for this transcript
	def mappability(self):
		minter = 1
		mintra = 1

		if self.num_inter > 0:
			# we have some inter-gene hits
			hlist, hvals, hrats = self.hits_to_lists(type="inter")
			minter = 1-np.percentile(hrats, 95)

		if self.num_intra > 0:
			# we have some intra-gene hits
			hlist, hvals, hrats = self.hits_to_lists(type="intra")
			mintra = 1-np.percentile(hrats, 95)


		return [mintra, minter]

	def pheader(self):

		lout = [
			"transcript_id", "gene_id", "gene_name", 
			"intra_mapp", "intra_min_mapp", "inter_mapp", "inter_min_mapp",
			"length", "num_kmers", 
			"intra_transcripts_hit", "inter_transcripts_hit", 
			"intra_tid_list", "intra_tid_sharedK", 
			"inter_tid_list", "inter_tid_sharedK"
		]

		return("\t".join(lout))

	def tolist(self):

		intra_hlist = ""
		intra_hvals = ""
		inter_hlist = ""
		inter_hvals = ""
		intra_min = 1
		inter_min = 1

		if self.num_intra > 0:
			hlist, hvals, hrats = self.hits_to_lists(type="intra")
			intra_hlist = ";".join(hlist)
			intra_hvals = ";".join(map(str,hvals))
			intra_min = 1-max(hrats)

		if self.num_inter > 0:
			hlist, hvals, hrats = self.hits_to_lists(type="inter")
			inter_hlist = ";".join(hlist)
			inter_hvals = ";".join(map(str, hvals))
			inter_min = 1-max(hrats)

		# id
		lout = [self.id, self.gid, self.gname]

		# append inter and intra mappability
		mpp = self.mappability()
		
		lout += [ftos(mpp[0]), ftos(intra_min), ftos(mpp[1]), ftos(inter_min)]

		# transcript length
		lout.append(self.length)

		# number of kmers
		lout.append(self.num_k)

		# append number of hits for inter and intra
		lout.append(self.num_intra)
		lout.append(self.num_inter)

		# get lists of hits for intra
		if self.num_intra > 0:
			lout.append(intra_hlist)
			lout.append(intra_hvals)
		else:
			lout += ["u", "u"]

		if self.num_inter > 0:
			lout.append(inter_hlist)
			lout.append(inter_hvals)
		else:
			lout += ["u", "u"]

		return lout

def ftos(x):
	return("{:0.4f}".format(x))

def quantile(v, q):
	n = len(v)
	vhat = sorted(v)

#
# write results to file
def print_results(dtid, fname):

	fout = open(fname, "w")

	ktid = sorted(dtid.keys())

	fout.write(dtid[ktid[0]].pheader() + "\n")

	for tid in ktid:

		fout.write("\t".join(map(str, dtid[tid].tolist())) + "\n")

	fout.close()

	return 0


def parse_hits(ff, tid, dtid, dfai):

	#
	# setup entry for this transcript. we'll keep the transcript's own kmer
	# count and then shared counts with other transcripts. those will be kept
	# in a dict
	tout = Tobj(tid, length=dfai[tid])

	#
	# open hits
	fin = open(ff, "r")
	for szl in fin:
		if szl[0] != ">":
			continue

		aln = szl.strip().split("\t")
		for i in range(1, len(aln)):
			# extract target tid
			ttid = (aln[i].split("="))[0]
			if ttid==tid:
				tout.kpp()
			else:
				tout.add_inter_hit(ttid)

	fin.close()

	return tout


def parse_ii_hits(ff, tid, dtid, tid2gid, tid2gname, dfai):

	#
	# setup entry for this transcript. we'll keep the transcript's own kmer
	# count and then shared counts with other transcripts. those will be kept
	# in a dict
	tout = Tobj(tid, gid=tid2gid[tid], gname=tid2gname[tid], length=dfai[tid])

	#
	# open hits
	fin = open(ff, "r")
	for szl in fin:
		if szl[0] != ">":
			continue

		aln = szl.strip().split("\t")
		for i in range(1, len(aln)):
			# extract target tid
			ttid = (aln[i].split("="))[0]
			tgid = tid2gid[ttid]

			if ttid==tid:
				tout.kpp()
			else:
				if tgid == tout.gid:
					tout.add_intra_hit(ttid)
				else:
					tout.add_inter_hit(ttid)

	fin.close()

	return tout

def parse_info(fname):

	tid2gid = {}
	tid2gname = {}

	# open the file
	fin = open(fname, "r")
	for szl in fin:
		if re.search("transcript_id", szl):
			# header line, skip it
			continue

		aln = szl.strip().split("\t")
		tid = aln[3]

		tid2gid[tid] = aln[2]
		tid2gname[tid] = aln[4]

	fin.close()

	return tid2gid, tid2gname


def parse_fai(fai):
	dfai = {}
	fin = open(fai, "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		dfai[aln[0]] = int(aln[1])

	fin.close()
	return dfai

#
# check for fai index for the supplied fasta
def fai_exists(fa):
	return(isfile("{}.fai".format(fa)))

#
# if the supplied fasta file doesn't have an fai index then this function
# calls samtools to build it
def build_fai(fa):
	cmd = "samtools faidx {}".format(fa)
	runcmd(cmd)
	return 0

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
	sys.stderr.write("[kmer-transcript-mappability] " + sz + "\n")



#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")

# -- required args
parser.add_argument('fasta', type=str, help="Input FASTA")

# -- options
parser.add_argument('-i', '--info', default=None, action="store", 
	help="so-called '.info' file generated from a GTF. enables the program to differentiate between inter-gene and intra-gene mappabilities")
parser.add_argument('-t', '--threads', default=1, action="store", 
	help="Number of threads for seal [1]")
parser.add_argument('-k', '--kmer-length', default=31, type=int, action="store", 
	help="Kmer length for matching. Must be an odd number [31]")
parser.add_argument('-v', action="store_const", const=True, default=False, 
	help="Print more messages to stderr while running")

# -- parse command line
args = parser.parse_args()


if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

