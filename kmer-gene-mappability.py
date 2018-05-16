#!/usr/bin/python
#==============================================================================
# kmer-gene-mappability.py
#
# Shawn Driscoll
# 20170315
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Extracts kmers grouped by gene to determin each gene's mappability with
# respect to the entire transcriptome.  
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
	#
	# variables
	szfai = "{}.fai".format(args.fasta)
	rres = None
	tid2gid = {}
	tid2gname = {}
	gname2tid = {}
	gid2tid = {}
	has_info = False
	ntid = 0
	total_proc = 0
	outfile = ""
	
	if not isfile(args.fasta):
		message("Error: input FASTA file does not exist")
		return 1

	#
	# check for fasta index
	if not fai_exists(args.fasta):
		# need to build on
		message("building FAI index")
		build_fai(args.fasta)

	#
	# parse the fai. this is used to make sure genes have 
	# at least one transcript that is at least as long as the 
	# 'kmer_length' setting
	message("loading FAI index {}".format(szfai))
	dfai = parse_fai(szfai)
	
	##
	# create output file name
	if args.o is None:
		args.o = drop_file_ext(args.fasta) + ".gmapp"
	else:
		args.o = args.o + ".gmapp"
	
	
	#
	# load annotation
	message("loading {}".format(args.info))
	tid2gid, tid2gname, gname2tid, gid2tid = parse_info(args.info)
	has_info = True

	#
	# setup the Queue

	message("starting {} worker threads for kmer counting".format(args.threads))

	if args.gene_name:
		# process gene names
		kq = KmerQueue(args.threads, args.kmer_length, args.fasta, gname2tid)
		ttable = gname2tid
		titable = tid2gname
	else:
		kq = KmerQueue(args.threads, args.kmer_length, args.fasta, gid2tid)
		ttable = gid2tid
		titable = tid2gid

	#
	# get all ids for the following loop
	all_ids = sorted(ttable.keys())
	n = len(all_ids)
	# keep track of genes that were skipped. this list resets
	# at each iteration of the main loop
	skipped_too_short = []
	# two counters used to step through the full list in 
	# windows equal to args.m
	i0 = 0
	i = 0

	with open(args.o, "w") as fout:

		# main loop. exit condition is tested within
		while True:
	
			#
			# inner loop to grab 'args.m' genes. gene names are pushed 
			# to the queue which sends them out to the worker threads
	
			i = 0
			message("fetching {} genes and counting kmers...".format(min([args.m, n-i0])))
	
			while i < args.m:
	
				gid = all_ids[i0+i]
				tset = ttable[gid]
	
				# confirm we have some features that are long enough to get kmer'd
				num_short = 0
				for tid in tset:
					if dfai[tid] < args.kmer_length:
						# short sequence to be just printed back out
						# skipped_too_short.append(tid)
						num_short += 1
	
				if num_short < len(tset):
					# good process this gene
					kq.q.put(gid)
	
				elif num_short == len(tset):
					# nope, skipping the entire gene
					skipped_too_short.append(gid)
	
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
				cmd = "seal.sh in=tmp.mer.fa ref={} k={} threads={} rcomp=f mm=f hdist=0 out=tmp.hits.fa rename=t ambig=all tuc=t trd=t overwrite=t -Xmx{} 2>seal.log".format(args.fasta, args.kmer_length, args.seal_threads, args.x)
				runcmd(cmd, args.v)
				# parse hits
	
				rres = parse_hits("tmp.hits.fa", tid2gid, tid2gname, ttable, titable)
	
				runcmd("rm tmp.hits.fa", False)
	
			elif len(skipped_too_short) > 0:
				# we have some that are too short so we'll insert some blanks
				message("found {} sequences shorter than kmer length")
				for gid in skipped_too_short:
					tset = ttable[gid]
					rres[gid] = Tobj(gid, tset=ttable[gid])
					for tid in tset:
						rres[gid].gene_name.update([tid2gname[tid]])
						rres[gid].gene_id.update([tid2gid[tid]])
	
			if i0==0:
				# print header
				fout.write(rres[rres.keys()[0]].pheader() + "\n")
	
			for tt in rres.keys():
				fout.write(rres[tt].__str__() + "\n")
	
			i0 += args.m
			skipped_too_short = []
	
			message("processed {} of {} genes. {:0.1f}% complete".format(total_proc, n, total_proc*100.0/n))
	
			if i0 >= n:
				break

	if isfile("tmp.mer.fa"):
		runcmd("rm tmp.mer.fa")

	return 0


#==============================================================================
# defs
#==============================================================================

def drop_file_ext(fname):
	tmp = fname.split(".")
	return ".".join(tmp[0:(len(tmp)-1)])

#
# this class provides a multi-thread queue for processing sequences into mers.
# a single FASTA is written per thread. the names of the FASTA files will be
# <thread name>.fa and when done you can use the values in the 'self.names' 
# attribute to check them for content. 
class KmerQueue(object):
	def __init__(self, n, k, ref, ttable):
		
		self.n = n
		self.k = k
		self.ref = ref
		self.q = Queue()
		# table to translate gene names/ids to transcripts
		self.ttable = ttable 
		# set to hold names of the threads
		self.names = set()
		
		# start all threads
		for i in range(n):
			t = threading.Thread(target=self.worker)
			t.daemon = True
			t.start()

	def worker(self):
		name = threading.currentThread().getName()
		self.names.update([name])

		while True:
			item = self.q.get()
			# get transcript names 
			tset = " ".join(self.ttable[item])
			cmd = "samtools faidx {} {} | ".format(self.ref, tset)
			cmd += "kmercountexact.sh threads=1 in=stdin.fa k={} rcomp=f tuc=t trd=t overwrite=t out=stdout.fa 2>/dev/null | ".format(self.k)
			cmd = cmd + FASTA_NAME_PATH + " '{};' - 2>/dev/null >> {}.fa".format(item, name)
			runcmd(cmd, verbose=False)

			self.q.task_done()


class Tobj(object):
	def __init__(self, gid, tset):

		self.gid = gid
		self.tset = tset
		self.num_k = 0
		self.inter_hits = {}		
		self.num_inter = 0
		self.gene_id = set()
		self.gene_name = set()


	def __str__(self):
		lout = self.tolist()
		return "\t".join(map(str, lout))

	def add_inter_hit(self, tid):
		if tid not in self.inter_hits:
			self.inter_hits[tid] = 0
			self.num_inter += 1

		self.inter_hits[tid] += 1
		return 0

	# increment the kmer counter
	def kpp(self):
		self.num_k += 1

	def hits_to_lists(self):

		foo = True
		hlist = []
		hvals = []
		hrats = []

		if self.num_inter == 0:
			foo = False
		hits = self.inter_hits

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

		if self.num_inter > 0:
			# we have some inter-gene hits
			hlist, hvals, hrats = self.hits_to_lists()
			minter = 1-np.percentile(hrats, 95)

		return minter

	def pheader(self):

		lout = [
			"group_id", "gene_id", "gene_name", "transcript_id",  
			"mapp", "min_mapp",
			"num_kmers", 
			"num_genes_hit", 
			"gene_list", "gene_sharedK"
		]

		return("\t".join(lout))

	def tolist(self):

		inter_hlist = ""
		inter_hvals = ""
		inter_min = 1

		if self.num_inter > 0:
			hlist, hvals, hrats = self.hits_to_lists()
			inter_hlist = ";".join(hlist)
			inter_hvals = ";".join(map(str, hvals))
			inter_min = 1-max(hrats)

		# id
		lout = [self.gid, ";".join(sorted(list(self.gene_id))), 
			";".join(sorted(list(self.gene_name))), ";".join(self.tset)]

		# append inter and intra mappability
		mpp = self.mappability()
		
		lout += [ftos(mpp), ftos(inter_min)]

		# number of kmers
		lout.append(self.num_k)

		# append number of hits for inter and intra
		lout.append(self.num_inter)

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
# this function parses the seal output to collect all of the kmer counts 
# per feature that was processed. those are returned as a dict of Tobj
# objects. 
def parse_hits(ff, tid2gid, tid2gname, ta, tb):
	# ta translates gene names to transcript sets
	# tb translates transcript ids to gene name

	dd = {}

	#
	# open hits
	fin = open(ff, "r")
	for szl in fin:
		if szl[0] != ">":
			continue

		aln = szl.strip().split("\t")

		# 
		# gid is in the first element
		tmp = aln[0].split(";")
		gid = tmp[0][1:]
		
		if gid not in dd:
			# create new Tobj for this gene
			tset = ta[gid]
			dd[gid] = Tobj(gid, tset=ta[gid])
			# populate this thing's gene names and gene ids
			for tid in tset:
				dd[gid].gene_name.update([tid2gname[tid]])
				dd[gid].gene_id.update([tid2gid[tid]])

		gid_hit = set()
		for i in range(1, len(aln)):
			# extract target tid and translate it
			ttid = (aln[i].split("="))[0]
			tgid = tb[ttid]

			if tgid==gid:
				# hit to self, increment kmer count
				dd[gid].kpp()
			else:
				# hit to other. increment shared kmer count for that feature
				if tgid not in gid_hit:
					# per kmer (line from file) we can only count a hit to a gene id once
					# because a kmer may be shared between many transcripts of a gene
					dd[gid].add_inter_hit(tgid)
					gid_hit.update([tgid])

	fin.close()

	return dd

#
# parse_info
# This function loads the gene/transcript annotation into a 
# few dicts.
def parse_info(fname):

	# four dicts for holding info. 
	tid2gid = {}		# transcript to gene id
	tid2gname = {}		# transcript to gene name
	gid2tid = {}		# gene id to transcript set
	gname2tid = {}		# gene name to transcript set

	# open the file
	fin = open(fname, "r")
	for szl in fin:
		
		if re.search("transcript_id", szl):
			# header line, skip it
			continue

		# split the line up
		aln = szl.strip().split("\t")
		# gather info
		tid = aln[3]
		gid = aln[2]
		# make sure there aren't any spaces in the gene name
		gname = re.sub(" ", "_", aln[4])

		tid2gid[tid] = gid
		tid2gname[tid] = gname

		# add gene id to transcript set entry
		if gid not in gid2tid:
			gid2tid[gid] = []

		gid2tid[gid].append(tid)

		# add gene name to transcript set entry
		if gname not in gname2tid:
			gname2tid[gname] = []

		gname2tid[gname].append(tid)

	fin.close()

	return tid2gid, tid2gname, gname2tid, gid2tid


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
	sys.stderr.write("[kmer-gene-mappability] " + sz + "\n")



#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Generates distinct kmers by gene (pooling multiple transcripts if they exist) and then quantifies with which other genes the gene shares its kmers with.")

# -- required args
parser.add_argument('fasta', type=str, help="FASTA to analyze")
parser.add_argument('info', type=str, 
	help="so-called '.info' file generated from a GTF.")

# -- options
parser.add_argument('-o', type=str, action="store", default=None, 
	help="Output file name stub (gmapp) is appended automatically. Default is to use the input fasta file name.")
parser.add_argument('-N', '--gene-name', action="store_const", const=True, default=False, 
	help="Process by gene name instead of gene id")
parser.add_argument('-t', '--threads', default=2, type=int, action="store", 
	help="Number of threads for generating kmers [2]")
parser.add_argument('-p', '--seal-threads', default=1, type=int, action="store", 
	help="Number of threads for Seal [1]")
parser.add_argument('-k', '--kmer-length', default=31, type=int, action="store", 
	help="Kmer length for matching. Must be an odd number [31]")
parser.add_argument('-v', action="store_const", const=True, default=False, 
	help="Print more messages to stderr while running")
parser.add_argument('-m', default=100, type=int, action="store", 
	help="Number of genes to process together [100]")
parser.add_argument('-x', type=str, default="4g", 
	help="Enter custom value for seal for memory allocation (i.e. for -Xmx10g enter 10g [4g]")

# -- parse command line
args = parser.parse_args()


if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

