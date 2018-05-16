#!/usr/bin/python
#==============================================================================
# scrna-samhits.py
#
# Shawn Driscoll
# 20171208
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This code parses transcriptome alignments, expecting multiple 
# mappings per read, for a single cell and produces TCC type counts
# like kallisto pseudo. This is a UMI counter though it could be modified to 
# be a read counter. The UMI barcode must be appended to the read name
# with a ":" as a delimeter.
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime, time
import subprocess as sp
from os import system, unlink
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock
import hashlib
import igraph as ig

# from igraph import *
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

	# variables
	err = False
	sam = False
	read_orientation = 1
	if args.unstranded:
		read_orientation = 0

	##
	## confirm all files exist. check for target file.
	if not isfile(args.fin):
		error_message("Input file {} does not exist".format(args.fin))
		err = True

	if not isfile(args.ref):
		error_message("Input file {} does not exist".format(args.ref))
		err = True

	if err:
		return 1

	##
	## all good?
	##

	##
	## check file type
	if re.search("\.sam$", args.fin):
		sam = True

	##
	## load reference. refFlat format is super easy because the gene name 
	## and transcript ids are the first two columns.
	##
	if not args.quiet:
		message("Loading reference {}".format(args.ref))
	dref = {}
	dgenes = {}
	t0 = time()
	with open(args.ref, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			aln[0] = re.sub("\s", "_", aln[0])
			dref[aln[1]] = aln[0]
			if aln[0] not in dgenes:
				dgenes[aln[0]] = 0
	
	if not args.quiet:
		sys.stderr.write("{} sec\n".format(time()-t0))
	
	##
	## deal with input file. if it is bam then we have to write it out as sam
	## temporarily. we need sam format because we are gonna index the reads
	## and we'll need to have seek access.
	##
	
	tmpfile = hashlib.md5(args.fin).hexdigest()
	if not sam:
		if not args.quiet:
			message("Converting BAM to SAM to allow indexing")
		t0 = time()
		p1 = sp.Popen("samtools view -h -o {}.sam {}".format(tmpfile, args.fin).split())
		p1.wait()
		if not args.quiet:
			sys.stderr.write("{} sec\n".format(time()-t0))
		alignments = "{}.sam".format(tmpfile)
	else:
		alignments = args.fin
	
	t0 = time()
	if not args.quiet:
		message("Indexing reads...")
	read_index = build_index(alignments, read_orientation, args.q)
	if not args.quiet:
		sys.stderr.write("{} sec\n".format(time()-t0))

	##
	## build the EC's
	t0 = time()
	if not args.quiet:
		message("Building equivalence class table")
	ec_table = {}
	
	# loop through read index and make additional EC's for reads that hit multiple targets
	fin = open(alignments, "r")
	
	# this dict is key'd by unique UMIs. a single UMI may be associated with several ecid's
	# and only by disambiguating this can we uniquly assign the UMI to an ecid
	umi_table = defaultdict(list)
	
	for rid in read_index.keys():
		ecid = None
		
		if len(read_index[rid]) > 1:
			# multi-mapper. get the targets
			tidset = set()
			for offset in read_index[rid]:
				fin.seek(offset)
				aln = fin.readline().strip().split("\t")
				tid = dref[get_aligned_target(aln, 0)]
				tidset.add(tid)
			
			tidset = sorted(list(tidset)) 
			ecid = ":".join(tidset)
		
		else:
			# get this read's target
			fin.seek(read_index[rid][0])
			aln = fin.readline().strip().split("\t")
			ecid = dref[get_aligned_target(aln, 0)]
		
		if ecid is not None:
			# we have a target
			# get the UMI from the parsed read

			if ecid not in ec_table:
				ec_table[ecid] = []

			umi = parse_umi(aln[0])
			ec_table[ecid].append(umi)
			
			umi_table[umi].append(ecid)
	
	if not args.quiet:
		sys.stderr.write("{} sec\n".format(time()-t0))

	fin.close()

	final_counts = defaultdict(float)

	# initalize with all of the annotated gene names
	for gname in dgenes.keys():
		final_counts[gname] = 0

	for umi in umi_table.keys():
		ecset, eccount = count_umi(umi_table[umi])

		# check total count. we only want to keep UMI that have >= args.min_umi_freq
		# occurences
		if sum(eccount) < args.min_umi_freq:
			continue
			
		# get final list of targets for this UMI
		targets = cluster_ecid(ecset)
		# since we are doing this in terms of genes it's possible to get the same
		# UMI at more than one gene so we'll count the whole UMI to each assigned feature.
		w = 1.0
		for ecid in targets:
			final_counts[ecid] += w
	
	##
	## right here we can attempt to resolve hits to multiple features iterativly like 
	## I did in 'annotate-kallisto-pseudo.py'
	##
	
	# initalize weights and separate out the ids that pertain to only a single 
	# gene
	ecid_multi_weights = {}
	single_counts = {}
	for ecid in final_counts.keys():
		tmp = ecid.split(":")
		if len(tmp) > 1:
			n = len(tmp)
			ecid_multi_weights[ecid] = [1.0/n for i in range(n)]
		else:
			single_counts[ecid] = final_counts[ecid]
	
	num_multi = len(ecid_multi_weights.keys())
	if not args.quiet:
		sys.stderr.write("Found {} multi-gene ids\n".format(len(ecid_multi_weights.keys())))
	
	max_iter = 100
	iter = 0
	tol = 1e-6
	
	# used to calculate weights at each iteration
	hits_dev0 = copy_counts(single_counts)
	# used to update counts
	hits_dev = copy_counts(single_counts)
	
	sse0 = 0
	sse = 0
	
	if num_multi > 0:
		while True:
			iter += 1
			sse = 0
			
			for ecid in ecid_multi_weights.keys():
				w0 = ecid_multi_weights[ecid]
				gset = ecid.split(":")
				# get gene counts
				w = sumnorm([hits_dev0[g] for g in gset])
				# get mean squared error
				sse += sum([(w0[i]-w[i])**2 for i in range(len(w))])*1.0/len(w)
				# update weight
				ecid_multi_weights[ecid] = list(w)
				
				if sum(w) > 0:
					# print ecid, w0, w
					# assign hits
					for i in range(len(w)):
						hits_dev[gset[i]] += w[i]*final_counts[ecid]
				
			# update hits_dev0 and hits_dev
			hits_dev0 = copy_counts(hits_dev)
			hits_dev = copy_counts(single_counts)
			
			sse = sse*1.0/num_multi
			
			if not args.quiet:
				sys.stderr.write("Iteration: {}; mse: {}\n".format(iter, sse))
	
			if iter >= max_iter:
				if not args.quiet:
					sys.stderr.write("did not converge\n")
				break
			
			if abs(sse0-sse) < tol:
				if not args.quiet:
					sys.stderr.write("converged\n")
				break
			
			sse0 = sse
	
	# print total reads
	print "#total_reads={}".format(len(read_index.keys()))
	for ecid in sorted(hits_dev0.keys()):
		print "\t".join(map(str, [ecid, hits_dev0[ecid]]))

	return 0

def sumnorm(v):
	ss = sum(v)
	if ss > 0:
		vhat = [a*1.0/ss for a in v]
	else:
		vhat = list(v)
		
	return(vhat)

def copy_counts(v):
	x = {} 
	for id in v.keys():
		x[id] = v[id]
	return x

#
# given a list of distinct ecid's this clusters them and creates
# a final list of distict cluster ecids. for each cluster of ecids only
# the common transcript ids between them are kept. this is an attempt
# to narrow down, per umi, a distinct source transcript.
def cluster_ecid(ecset):
	
	ehat = []
	efinal = []
	n = len(ecset)
	
	if n == 1:
		return(ecset)
	
	# explode the ecids
	for ecid in ecset:
		ehat.append(ecid.split(":"))
	
	pairs = []
	
	# cluster
	
	for i in range(n-1):
		for j in range(i+1, n):
			if list_int(ehat[i], ehat[j]):
				pairs.append((i,j))
	
	if len(pairs) > 0:
		# need to cluster so use a graph to make this easier
		g = ig.Graph()
		g.add_vertices(n)
		g.add_edges(pairs)
		g_clust = g.clusters()
		bundles = g_clust.membership
		# bundles is a clustering membership vector
		num_clust = max(bundles)+1
		tmp = [[] for i in range(num_clust)]
		for i in range(n):
			tmp[bundles[i]] += ehat[i]
		
		# reduce the tid lists in each bundle down to only those that 
		# occure the maximum number of times
		for i in range(num_clust):
			tid, tid_c = count_umi(tmp[i])
			nmax = max(tid_c)
			tmp2 = []
			for j in range(len(tid)):
				if tid_c[j]==nmax:
					# keep this tid
					tmp2.append(tid[j])
			
			efinal.append(":".join(sorted(tmp2)))
	
	else:
		# no combinations are possible so we have to keep the original list
		efinal = list(ecset)
	
	return efinal
		
def list_int(a, b):
	rres = set(a).intersection(set(b))
	return len(rres) > 0

##
## given a list of umi this produces a list of unique umi with counts
def count_umi(lumi):
	
	dumi = {}
	umi = []
	umin = []
	
	for u in lumi:
		if u not in dumi:
			dumi[u] = 0
			umi.append(u)
			umin.append(0)
		
		dumi[u] += 1
	
	# transfer counts to the list
	for i in range(len(umi)):
		umin[i] = dumi[umi[i]]
	
	return umi, umin
		
	
##
# build_index
# this function indexes aligned reads. You can pass it a strand orientation so 
# only reads with that orientation are indexed:
# strand == 0: unstranded
# strand == 1: forward stranded (default)
# strand == -1: reverse stranded
def build_index(alignments, strand, min_mapq):
	
	offset = 0
	read_index = defaultdict(list)
	
	with open(alignments, "r") as fin:
		
		for szl in fin:
			
			if szl[0] != "@":
				# not a header line
				aln = szl.strip().split("\t")
				if int(aln[1]) & 0x4:
					# unaligned, skip it
					offset += len(szl)
					continue
				elif len(aln) < 10:
					# probably not a real alignment
					offset += len(szl)
					continue
				elif int(aln[4]) < min_mapq:
					# mapq too low
					offset += len(szl)
					continue
				
				elif strand != 0:
					read_orientation = 1
					if int(aln[1]) & 0x10:
						# alignment is reverse strand
						read_orientation = -1
					if read_orientation == strand:
						# keep it
						read_index[aln[0]].append(offset)
						offset += len(szl)
						continue
					else:
						# nope
						offset += len(szl)
						continue
				else:
					# keep it
					read_index[aln[0]].append(offset)
			
			# we get here if input line was a header line or if
			# kept the read.					
			offset += len(szl)
	
	return read_index



##
# this function returns the cell barcode and umi barcode from a read name string
def parse_umi(sz):
	tmp = sz.split(":")
	return tmp[-1]

##
# returns the target of an alignment depending on strand expectation. if the
# strand is expected to be same-strand then strand_type==1.  opposite is -1 and
# strand_type==0 means ignore strand.
def get_aligned_target(aln, strand_type):
	
	tid = None
	
	if (int(aln[1]) & 0x4):
		# unaligned
		return tid
	
	orientation = -1
	if (int(aln[1]) & 0x10) == 0:
		orientation = 1

	if strand_type != 0:
		if strand_type == orientation:
			tid = aln[2]
	
	else:
		tid = aln[2]
	
	return tid	

##
# runs bbduk to filter reads against a reference keeping any read with 
# at least a single kmer match
def bbduk_filter(infile, outfile, ref, k):

	cmd = "bbduk.sh in={} ref={} k={} hdist=0 mm=f overwrite=t".format(infile, ref, k)
	cmd += " outm={}".format(outfile)
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	
	return 0

##
# runs bbmap to filter reads keeping only those that pass some minimum 
# alignment identity
def bbmap_filter(infile, outfile, ref, identity):
	
	cmd = "bbmap.sh in={} path={} idfilter={} outm={} overwrite=t".format(infile, ref, identity, outfile)
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	
	return 0

##
# filters reads based on umi frequency. all distinct umis are counted and then 
# reads are written to 'outfile' that have UMIs which occur at least the 
# specified minimum number of times.
def umi_filter(infile, outfile, min_count, fmt, gz):

	# if data is gzipped then first we will undo that
	if gz:
		cmd = "gunzip {}".format(infile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)
		infile = re.sub("\.gz", "", infile)

	if fmt == "fasta":
		rres = umi_filter_fasta(infile, outfile, min_count)
	else:
		rres = umi_filter_fastq(infile, outfile, min_count)

	if gz:
		cmd = "gzip {}".format(infile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)
		cmd = "gzip {}".format(outfile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)

	return rres


def umi_filter_fasta(infile, outfile, min_count):

	idx = 0
	dumi = defaultdict(int)
	keep_read = False
	keep_count = 0
	drop_count = 0

	##
	# first count the UMIs
	with open(infile, "r") as fin:

		for szl in fin:
			if szl[0] == ">":
				# read name, get barcodes
				bc = parse_barcodes(szl.strip())
				dumi[bc[1]] += 1

	##
	# read back through and print out the reads passing the filter. also write
	# the umi file.
	szout = ""
	umiout = ""
	with open(infile, "r") as fin:

		for szl in fin:

			if szl[0] == ">":
				keep_read = False
				bc = parse_barcodes(szl.strip())
				if dumi[bc[1]] >= min_count:
					umiout += "{}\n".format(bc[1])
					keep_read = True
					keep_count += 1
				else:
					drop_count += 1

			if keep_read:
				szout += szl

	with open(outfile, "w") as fout:
		fout.write(szout)
	
	tmp = outfile.split(".")
	tmp[-1] = "umi"
	umifile = ".".join(tmp)
	with open(umifile, "w") as fout:
		fout.write(umiout)
	

	return [keep_count, drop_count]


def umi_filter_fastq(infile, outfile, min_count):

	idx = 0
	dumi = defaultdict(int)
	keep_read = False
	keep_count = 0
	drop_count = 0

	##
	# first count the UMIs
	with open(infile, "r") as fin:

		for szl in fin:
			idx += 1
			if idx == 1:
				# read name, get barcodes
				bc = parse_barcodes(szl.strip())
				dumi[bc[1]] += 1

			if idx == 4:
				idx = 0

	##
	# read back through and print out the reads passing the filter
	szout = ""
	umiout = ""
	idx = 0
	with open(infile, "r") as fin:

		for szl in fin:
			idx += 1
			if idx == 1:
				keep_read = False
				bc = parse_barcodes(szl.strip())
				if dumi[bc[1]] >= min_count:
					umiout += "{}\n".format(bc[1])
					keep_read = True
					keep_count += 1
				else:
					drop_count += 1

			if keep_read:
				szout += szl

			if idx == 4:
				idx = 0

	with open(outfile, "w") as fout:
		fout.write(szout)

	tmp = outfile.split(".")
	tmp[-1] = "umi"
	umifile = ".".join(tmp)
	with open(umifile, "w") as fout:
		fout.write(umiout)

	return [keep_count, drop_count]


def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0


def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def error_message(sz):
	sys.stderr.write("Error: {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def print_exception():
	exc_type, exc_value, exc_traceback = sys.exc_info()
	print "*** print_tb:"
	traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
	print "*** print_exception:"
	traceback.print_exception(exc_type, exc_value, exc_traceback,
	                          limit=2, file=sys.stdout)
	print "*** print_exc:"
	traceback.print_exc()
	print "*** format_exc, first and last line:"
	formatted_lines = traceback.format_exc().splitlines()
	print formatted_lines[0]
	print formatted_lines[-1]
	print "*** format_exception:"
	print repr(traceback.format_exception(exc_type, exc_value,
	                                      exc_traceback))
	print "*** extract_tb:"
	print repr(traceback.extract_tb(exc_traceback))
	print "*** format_tb:"
	print repr(traceback.format_tb(exc_traceback))
	print "*** tb_lineno:", exc_traceback.tb_lineno	

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('ref', type=str, help="refFlat annotation to relate transcript ids to gene names")
parser.add_argument('fin', type=str, help="Genome aligned reads in sam/bam format")
parser.add_argument('--min-umi-freq', type=int, default=1, 
	help="Minimum frequency of a UMI within a cell for output to file.")
parser.add_argument("-q", type=int, default=0, 
	help="Minimum MAPQ for accepting alignments.")
parser.add_argument("--quiet", action="store_const", const=True, default=False, 
	help="Quiet mode. No messages are printed out.")
parser.add_argument("--unstranded", action="store_const", const=True, default=False, 
	help="As of this code's writing we run the 10x libraries stranded. Second mate is the read we align and that one would always be in forward orientation. Set this flag if that is not true.")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

