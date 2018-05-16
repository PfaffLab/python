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
from Queue import Queue
import threading

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
FASTA_NAME_PATH = "{}/coding/perl/fasta-name-prefix.pl".format(HOME)

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
	total_proc = 0

	#
	# check for fasta index
	if not fai_exists(args.fasta):
		# need to build on
		message("building FAI index")
		build_fai(args.fasta)

	#
	# parse the fai
	message("loading FAI index {}".format(szfai))
	dfai = parse_fai(szfai)
	
	#
	# if info is supplied load it
	if args.info is not None:
		message("loading {}".format(args.info))
		tid2gid, tid2gname = parse_info(args.info)
		has_info = True

	#
	# get all tid keys
	#

	all_tid = sorted(dfai.keys())
	n = len(all_tid)
	skipped_too_short = []
	i0 = 0
	i = 0

	message("starting {} worker threads for kmer counting".format(args.threads))
	kq = KmerQueue(args.threads, args.kmer_length, args.fasta)

	while True:

		#
		# fetch sequences, kmer and alter their sequence names to include the transcript id
		i = 0
		message("fetching {} sequences and counting kmers...".format(min([args.m, n-i0])))

		while i < args.m:

			tid = all_tid[i0+i]

			if dfai[tid] >= args.kmer_length:
				# add tid to the buffer
				kq.q.put(tid)
			else:
				# short sequence to be just printed back out
				skipped_too_short.append(tid)

			i += 1
			total_proc += 1

			if i0+i >= n:
				# force the loop to exit
				i = args.m

		# finish queue
		kq.q.join()
		
		# match the buffered kmers to the full index

		if len(skipped_too_short) < args.m:
			#
			# number of skipped is less than the total that were parsed out so we have stuff to match
			# and process
			#

			# merge the queue outputs
			runcmd("cat /dev/null > tmp.mer.fa", False)
			for tname in list(kq.names):
				if isfile("{}.fa".format(tname)):
					runcmd("cat {}.fa >> tmp.mer.fa".format(tname), False)
					runcmd("rm {}.fa".format(tname), False)

			message("running seal...")
			# seal is run without reverse complement matching and not allowed errors. all ambiguous hits are 
			# printed
			cmd = "seal.sh in=tmp.mer.fa ref={} k={} threads={} rcomp=f mm=f hdist=0 out=tmp.hits.fa rename=t ambig=all tuc=t trd=t overwrite=t 2>seal.log".format(args.fasta, args.kmer_length, args.seal_threads)
			runcmd(cmd, args.v)
			# parse hits

			message("parsing hits...")
			if has_info:
				rres = parse_ii_hits("tmp.hits.fa", tid2gid, tid2gname, dfai)
			else:
				rres = parse_hits("tmp.hits.fa", dfai)

			runcmd("rm tmp.hits.fa", False)

		elif len(skipped_too_short) > 0:
			# we have some that are too short so we'll insert some blanks
			message("found {} sequences shorter than kmer length")
			for tid in skipped_too_short:
				if has_info:
					rres[tid] = Tobj(tid, gid=tid2gid[tid], gname=tid2gname[tid], length=dfai[tid])
				else:
					rres[tid] = Tobj(tid, length=dfai[tid])

		if i0==0:
			# print header
			print rres[rres.keys()[0]].pheader()

		for tt in rres.keys():
			print rres[tt]

		i0 += args.m
		skipped_too_short = []

		message("processed {} of {} sequences. {:0.1f}% complete".format(total_proc, n, total_proc*100.0/n))

		if i0 >= n:
			break

	if isfile("tmp.mer.fa"):
		runcmd("rm tmp.mer.fa")

	return 0


#==============================================================================
# defs
#==============================================================================

#
# this class provides a multi-thread queue for processing sequences into mers.
# a single FASTA is written per thread. the names of the FASTA files will be
# <thread name>.fa and when done you can use the values in the 'self.names' 
# attribute to check them for content. 
class KmerQueue(object):
	def __init__(self, n, k, ref):
		
		self.n = n
		self.k = k
		self.ref = ref
		self.q = Queue()
		self.names = set()
		
		for i in range(n):

			t = threading.Thread(target=self.worker)
			t.daemon = True
			t.start()

	def worker(self):
		name = threading.currentThread().getName()
		self.names.update([name])

		while True:
			item = self.q.get()

			cmd = "samtools faidx {} {} | ".format(self.ref, item)
			cmd += "kmercountexact.sh threads=1 in=stdin.fa k={} rcomp=f tuc=t trd=t overwrite=t out=stdout.fa 2>/dev/null | ".format(self.k)
			cmd = cmd + FASTA_NAME_PATH + " '{};' - 2>/dev/null >> {}.fa".format(item, name)
			runcmd(cmd, verbose=False)

			self.q.task_done()


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


def parse_hits(ff, dfai):

	#
	# setup entry for this transcript. we'll keep the transcript's own kmer
	# count and then shared counts with other transcripts. those will be kept
	# in a dict
	dd = {}

	#
	# open hits
	fin = open(ff, "r")
	for szl in fin:
		if szl[0] != ">":
			continue

		aln = szl.strip().split("\t")

		# 
		# tid is in the first element
		tmp = aln[0].split(";")
		tid = tmp[0][1:]

		if tid not in dd:
			dd[tid] = Tobj(tid, length=dfai[tid])

		for i in range(1, len(aln)):
			# extract target tid
			ttid = (aln[i].split("="))[0]
			if ttid==tid:
				dd[tid].kpp()
			else:
				dd[tid].add_inter_hit(ttid)

	fin.close()

	return dd


def parse_ii_hits(ff, tid2gid, tid2gname, dfai):

	#
	# setup entry for this transcript. we'll keep the transcript's own kmer
	# count and then shared counts with other transcripts. those will be kept
	# in a dict
	# tout = Tobj(tid, gid=tid2gid[tid], gname=tid2gname[tid], length=dfai[tid])

	dd = {}

	#
	# open hits
	fin = open(ff, "r")
	for szl in fin:
		if szl[0] != ">":
			continue

		aln = szl.strip().split("\t")

		# 
		# tid is in the first element
		tmp = aln[0].split(";")
		tid = tmp[0][1:]

		if tid not in dd:
			dd[tid] = Tobj(tid, gid=tid2gid[tid], gname=tid2gname[tid], length=dfai[tid])

		for i in range(1, len(aln)):
			# extract target tid
			ttid = (aln[i].split("="))[0]
			tgid = tid2gid[ttid]

			if ttid==tid:
				dd[tid].kpp()
			else:
				if tgid == dd[tid].gid:
					dd[tid].add_intra_hit(ttid)
				else:
					dd[tid].add_inter_hit(ttid)

	fin.close()

	return dd

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
parser.add_argument('-t', '--threads', default=2, type=int, action="store", 
	help="Number of threads for generating kmers [2]")
parser.add_argument('-p', '--seal-threads', default=1, type=int, action="store", 
	help="Number of threads for Seal [1]")
parser.add_argument('-k', '--kmer-length', default=31, type=int, action="store", 
	help="Kmer length for matching. Must be an odd number [31]")
parser.add_argument('-v', action="store_const", const=True, default=False, 
	help="Print more messages to stderr while running")
parser.add_argument('-m', default=10, type=int, action="store", 
	help="Number of transcripts to process together [10]")

# -- parse command line
args = parser.parse_args()


if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

