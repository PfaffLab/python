#!/usr/bin/python
#==============================================================================
# scrna-seal-count.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Possible seal command: 
# bbrun seal -in <fastq_with_umi_in_readname> -ref <fasta> -k 31 -hdist 0 -mm f \
#	-rename -outm mu.fa -rcomp f -overwrite -trd -ambig all
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from os.path import isfile, expanduser, basename, dirname
from collections import defaultdict, Counter
from time import localtime, time
from multiprocessing import Pool
from os import system
from hashlib import md5
import numpy as np

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

np.random.seed(123)

HOME = expanduser("~")

# global table so that it may be accessed in a Pool map call function
UMI_TABLE = {}

# GTF fields
GTF_RNAME = 0
GTF_SOURCE = 1
GTF_FEATURE = 2
GTF_START = 3
GTF_END = 4
GTF_SCORE = 5
GTF_STRAND = 6
GTF_FRAME = 7
GTF_ATTRIBUTE = 8

#==============================================================================
# main
#==============================================================================

def mainff(args):

	sz = "CACATGCATC"

	print hamset1(sz)



def main(args):

	# variables
	
	#UMI_TABLE = {}
	umi_merged = {}

	total_umi = 0
	total_merged_umi = 0
	total_reads = 0
	corrected_umi = set()
	seal_file = ""
	remove_seal_file = False
	output_stub = ""
	output_path = ""

	gcounts1 = {}
	gcounts2 = {}

	gene_umis = {}

	#
	# Check if we need to run seal or just get straight to parsing hits
	#

	if(args.q is not None):
		# we're gonna run seal and the input file is a fastq file...yes?
		if(args.ref is None):
			# we can't do anything without a reference!
			sys.stderr.write("Error: if running Seal then you must specify '-r' to provide a reference")
			return(1)

		# good, run it
		seal_file = run_seal(args.q, args.ref, args.k, args.d)
		remove_seal_file = True

		output_path = dirname(args.q)
		output_stub = basename(args.q).split(".")[0]
		

	else:
		# not running seal so the input should be all we need to get started
		seal_file = args.a

		output_path = dirname(args.a)
		output_stub = basename(args.a).split(".")[0]


	# 
	# build the output file names
	#
	out_hits = output_stub + ".umi_counts"
	out_classes = output_stub + ".gene_classes"
	out_iterations = output_stub + ".iterations_mm"
	out_genes = output_stub + ".iterations_genes"
	if len(output_path) > 0:
		out_hits = output_path + "/" + out_hits
		out_classes = output_path + "/" + out_classes
		out_iterations = output_path + "/" + out_iterations
		out_genes = output_path + "/" + out_genes

	# check if these files exist

	if isfile(out_hits):
		if args.f:
			message("Output file {} exists and will be overwritten".format(out_hits))
		else:
			message("Output file {} exists!".format(out_hits))
			return(1)

	if isfile(out_classes):
		if args.f:
			message("Output file {} exists and will be overwritten".format(out_classes))
		else:
			message("Output file {} exists!".format(out_classes))
			return(1)

	#
	# now we can get on with the counting process
	#
	message("Loading {}".format(args.ref_flat))
	t0 = time()
	dannot, dgname2tid = parse_refflat(args.ref_flat)
	print time()-t0, "sec"

	message("Indexing seal output")
	line_offsets = []
	offset = 0
	t0 = time()
	with open(seal_file, "r") as fin:
		for szl in fin:
			if szl[0]==">":
				line_offsets.append(offset)
			
			offset += len(szl)

	print time()-t0, "sec"
	
	# total reads
	total_reads = len(line_offsets)
	message("Total reads: {}".format(total_reads))
	num_iter = 1
	if args.subsample > 0 or args.bootstrap:
		num_iter = args.iter
	
	# build gene name index
	gname_index = defaultdict(int)
	num_genes = 0
	for gname in sorted(dgname2tid.keys()):
		gname_index[gname] = num_genes
		num_genes += 1

	# final counts table for each iteration
	iter_counts = [[] for i in range(num_iter)]
		
#	message("Processing seal output")
	iter_index = 0
	
	for iter_index in range(num_iter):
		
		if args.bootstrap:
			# full depth bootstrap sample
			sset = np.random.choice(total_reads, total_reads)
		elif args.subsample > 0:
			if args.subsample < 1:
				tmp = math.floor((total_reads*1.0)*args.subsample)
			else:
				tmp = args.subsample
				tmp = min([total_reads, tmp])
			sset = np.random.choice(total_reads, int(tmp), replace=False)
		else:
			# all reads, no resampling
			sset = range(total_reads)
		
		UMI_TABLE.clear()
		gene_umis.clear()
		total_umi = 0
		total_sampled_reads = len(sset)
		total_mapped_reads = 0
		message("Parsing seal output")
		t0 = time()
	
		# 
		# this loop reads in the alignments. for each alignment associate the read's UMI
		# with all genes that it hit equally well. we build a dict 'UMI_TABLE' which is 
		# keyed by UMIs from the reads and holds a list of lists of gene names. each 
		# list of gene names is the list of genes that a read with that UMI hit equally 
		# well.
		with open(seal_file, "r") as fin:
	
			for offset_idx in sorted(sset):
				fin.seek(line_offsets[offset_idx])
				szl = fin.readline()

				# split it up
				ll = szl.strip().split("\t")
				if len(ll) < 2:
					# no alignment
					continue

				total_mapped_reads += 1
				
				# split up the read name
				rl = ll[0].split(":")
				umi = rl[-1]

				# split up the hits into tid/kmer_count pairs
				hits = [a.split("=") for a in ll[1:len(ll)]]
				# convert to a genes
				gidx = 0
				didx = {}
				gnames = []
				gkmer = []
				for k in hits:
					if k[0] in dannot:
						gname = dannot[k[0]]
						if gname not in didx:
							didx[gname] = gidx
							gkmer.append(float(k[1]))
							gnames.append(gname)
							gidx += 1
						else:
							gkmer[didx[gname]] = max([gkmer[didx[gname]], float(k[1])])
				
				if len(gkmer)==0:
					continue
				
				# get list of the best of the best
				try:
					best_k = max(gkmer)
				except:
					print ll
					sys.exit(1)
				best_idx = [i for i, j in enumerate(gkmer) if j==best_k]

				gbest = [gnames[i] for i in best_idx]
				if umi not in UMI_TABLE:
					UMI_TABLE[umi] = []

				UMI_TABLE[umi].append(gbest)

		print time()-t0, "sec"
		
		#
		# this loop narrows down each of the lists of lists of target genes per 
		# umi by joining overlapping gene sets. this preserves the possibility that
		# a UMI may have been used at more than one gene loci and keeps the counts
		# separate. the length of the list at each UMI_TABLE[umi] index tells us 
		# how many reads that UMI came from. when sublists of joined they are just 
		# concatenated. after merging we visit each umi and decide for each sublist
		# of target genes which genes to assign it to. we count the number of 
		# occurences of each gene name in the sublist and assign the UMI to those
		# equal to the highest count.  the umi is assigned to each gene and the 
		# read count given to each gene is the read count divided by  the 
		# number of genes we're assigning it to. 
		message("Generating candidate gene lists for each UMI")
		t0 = time()
		for umi in UMI_TABLE.keys():
			merged = True
			ll = UMI_TABLE[umi]
			while merged:
				merged = False
				for i in range(len(ll)-1):
					for j in range((i+1), len(ll)):
						iset = set(ll[i]).intersection(set(ll[j]))
						if len(iset) > 0:
							# join them
							ll[i] += ll[j]
							ll[j] = []
							merged = True

			# assign the UMI out to the genes
			for i in range(len(ll)):
				if len(ll[i])==0:
					continue
				
				genes = []
				counts = []
				for k, v in Counter(ll[i]).iteritems():
					genes.append(k)
					counts.append(v)
				
				# find best
				cmax = max(counts)
				cmax_idx = [i for i, j in enumerate(counts) if j==cmax]
				
				for k in cmax_idx:

					if genes[k] not in gene_umis:
						gene_umis[genes[k]] = { 'umi':[], 'count':[] }
					
					gene_umis[genes[k]]['umi'].append(umi)
					gene_umis[genes[k]]['count'].append(counts[k]*1.0/len(cmax_idx))
		
		print time()-t0, "sec"
		
		#
		# in this loop we merge UMIs assigned to a gene that are at 
		# most 1 base off from each other and satisfy a condition of one of them
		# having at least 2x the read count of the other. we also filter out 
		# UMIs with fewer than args.N reads AFTER merging. 

		message("Assigning raw and corrected UMI counts to genes")
		t0 = time()
		num_genes_detected = 0
		raw_umi_total = 0
		corrected_umi_total = 0
		# initialize the counts table for this iteration
		iter_counts[iter_index] = [[0,0] for i in range(num_genes)]
		for gname, d in gene_umis.iteritems():
			# deal with collapsing umi and creating final counts
			
			raw_count = sum([a >= args.N for a in d['count']])
			raw_umi_total += raw_count
			
			utemp = [list(a) for a in sorted(zip(d['umi'], d['count']), key=lambda x: x[1], reverse=True)]
			
			if len(utemp) > 1:
				# see if we can merge any of the umi
				merged = True
				while merged:
					merged = False
					for i in range(len(utemp)-1):
						if utemp[i] is None:
							continue
						for j in range((i+1), len(utemp)):
							if utemp[j] is None: 
								continue
							di = hamming_dist(utemp[i][0], utemp[j][0])
							if utemp[i][1] > (utemp[j][1]*2) and di < 2:
								# merge them
								merged = True
								utemp[i][1] += utemp[j][1]
								utemp[j] = None

				corrected_count = 0
				for i in range(len(utemp)):
					if utemp[i] is None:
						continue
					if utemp[i][1] >= args.N:
						corrected_count += 1

			else:
				
				corrected_count = raw_count
			
			corrected_umi_total += corrected_count
			
			# add counts for this gene
			iter_counts[iter_index][gname_index[gname]][0] = raw_count
			iter_counts[iter_index][gname_index[gname]][1] = corrected_count
			
			if corrected_count > 0:
				num_genes_detected += 1
			

		print time()-t0, "sec"

		message("> Total sampled reads:   {}".format(total_sampled_reads), show_date=False)
		message("> Total mapped reads:    {}".format(total_mapped_reads), show_date=False)
		message("> Total raw UMI:         {}".format(raw_umi_total), show_date=False)
		message("> Total corrected UMI:   {}".format(corrected_umi_total), show_date=False)
		message("> Total genes detected:  {}".format(num_genes_detected), show_date=False)

	#
	# finished reading the seal output file. if it was generated in this script call
	# then we can remove it now.
	if remove_seal_file:
		system("rm -f {}".format(seal_file))

	# calculate uncertainty if we had more than 1 iteration
	if num_iter > 1:
		message("Computing UMI count uncertainty")
		num_genes_detected = 0
		t0 = time()
		umi_umean = [0 for i in range(num_genes)]
		umi_uvar = [0 for i in range(num_genes)]
		umi_mean = [0 for i in range(num_genes)]
		umi_var = [0 for i in range(num_genes)]
		
		for i in range(num_genes):
			# collect all iteration counts for gene
			tmp = [0 for j in range(num_iter)]
			tmp_u = [0 for j in range(num_iter)]

			for iter_index in range(num_iter):
				
				tmp_u[iter_index] = iter_counts[iter_index][i][0]
				tmp[iter_index] = iter_counts[iter_index][i][1]
			
			umi_umean[i] = np.mean(tmp_u)
			umi_uvar[i] = np.var(tmp_u, ddof=1)
			umi_mean[i] = np.mean(tmp)
			umi_var[i] = np.var(tmp, ddof=1)
		
		print time()-t0, "sec"
		
		szout = "gene_name\traw_hits\tcorrected_hits\traw_hits_err\tcorrected_hits_err\n"
		for gname in sorted(dgname2tid.keys()):
			gidx = gname_index[gname]
			szout += gname
			szout += "\t{:0.4f}".format(umi_umean[gidx])
			szout += "\t{:0.4f}".format(umi_mean[gidx])
			szout += "\t{}".format(umi_uvar[gidx])
			szout += "\t{}".format(umi_var[gidx])
			szout += "\n"
	
			if umi_mean[gidx] > 0:
				num_genes_detected += 1

		with open(out_hits, "w") as fout:
			fout.write(szout)

		# I also want to export the individual iteration quantifications as a long table
		szout = "{}\t{}\t0\n".format(num_genes, num_iter)
		for iter_index in range(num_iter):
			for gname in sorted(dgname2tid.keys()):
				gidx = gname_index[gname]
				lout = [gidx, iter_index, iter_counts[iter_index][gidx][1]]
				szout += "\t".join(map(str, lout))
				szout += "\n"

		with open(out_iterations, "w") as fout:
			fout.write(szout)

		# write genes
		temp = []
		for gname in dgname2tid.keys():
			temp.append([gname, gname_index[gname]])
		temp.sort(key=lambda x: x[1])
		with open(out_genes, "w") as fout:
			for i in range(len(temp)):
				fout.write("\t".join(map(str, temp[i])))
				fout.write("\n")

		
		message("> Final genes detected:  {}".format(num_genes_detected), show_date=False)
		

	else:
		
		szout = "gene_name\traw_hits\tcorrected_hits\n"
		for gname in sorted(dgname2tid.keys()):
			gidx = gname_index[gname]
			szout += gname
			szout += "\t"
			szout += "{}".format("\t".join(map(str, iter_counts[0][gidx])))
			szout += "\n"
	
			if iter_counts[0][gidx][1] > 0:
				num_genes_detected += 1
	
		with open(out_hits, "w") as fout:
			fout.write(szout)
		

	return 0
		


#
# this function is called from map to blast through all UMIs to find 
# potential merges. Merges are 1-base off from the parent and the parent
# must have more reads supporting it than the child. Also they have to 
# have an intersecting list of target transcripts.
def pair_umi(bc):

	rres = []
	if bc is not None:
		neighbors = hamset1(bc)
		for neighbor in neighbors:
			if neighbor in UMI_TABLE:
				iset = set(UMI_TABLE[bc]['kcount'].keys()).intersection(UMI_TABLE[neighbor]['kcount'].keys())

				if len(iset) > 0:
					if UMI_TABLE[bc]['count'] > UMI_TABLE[neighbor]['count']:
						rres.append([bc, neighbor])

	return rres

def pair_umi0(bc):

	rres = []
	if bc is not None:
		neighbors = hamset1(bc)
		for neighbor in neighbors:
			if neighbor in UMI_TABLE:
				iset = UMI_TABLE[bc]['targets'].intersection(UMI_TABLE[neighbor]['targets'])

				if len(iset) > 0:
					if UMI_TABLE[bc]['count'] > UMI_TABLE[neighbor]['count']:
						rres.append([bc, neighbor])

	return rres

#
# this function takes the string sz and generates all 1-base hamming distance
# neighbors to it. 
def hamset1(sz):
	lets = list("ACTG")
	lsz = list(sz)

	cousins = []
	for i in range(len(lsz)):
		ch = lsz[i]
		for ch0 in lets:
			if ch0 != ch:
				ltmp = list(lsz)
				ltmp[i] = ch0
				cousins.append("".join(ltmp))

	return cousins

def hamming_dist(a, b):
	
	return sum([i != j for i, j in zip(a, b)])


def parse_refflat(f):
	tid2gname = {}
	gname2tid = defaultdict(list)

	with open(f, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			# tid is aln[1], gname is aln[0]
			tid2gname[aln[1]] = aln[0]
			gname2tid[aln[0]].append(aln[1])

	return tid2gname, gname2tid


# bbrun seal -in <fastq_with_umi_in_readname> -ref <fasta> -k 31 -hdist 0 -mm f -rename -outm mu.fa -rcomp f -overwrite -trd -ambig all

def run_seal(fastq, ref, k, hdist):
	target = md5(fastq).hexdigest() + ".fasta"
	cmd = "seal.sh in={} ref={} k={} hdist={} mm=f rename out={} outu={} rcomp=f overwrite trd ambig=all".format(fastq, ref, k, hdist, "seal_mapped.fasta", "seal_unmapped.fasta")
	message("CMD: {}".format(cmd))
	system(cmd)
	cmd = "cat seal_mapped.fasta seal_unmapped.fasta | shuffle.sh in=stdin.fa out={} overwrite=t".format(target)
	system(cmd)
	cmd = "rm -f seal_mapped.fasta seal_unmapped.fasta"
	system(cmd)
	return target

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0

def message(sz, show_date=True):
	if show_date:
		msg = "[{}] {}".format(time_string(), sz)
	else:
		msg = sz

	sys.stderr.write(sz + "\n")

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

parser = argparse.ArgumentParser(description="Count UMI hits to genes from Seal output. UMI barcodes are expected to be at the end of each read name separated by a single ':'.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('ref_flat', type=str, help="Annotation in refFlat format (use gtf2refFlat.py)")

finput = parser.add_mutually_exclusive_group(required=True)
finput.add_argument("-q", type=str, default=None, 
	help="Input FASTQ reads to process with Seal and then parse for UMI counts (requires -r)")
finput.add_argument("-a", type=str, default=None,
	help="Input FASTA file (with rename=t) from Seal. Count hits from this file directly")

parser.add_argument("-f", action="store_const", const=True, default=False, 
	help="Force overwrite of output files if they exist already")
parser.add_argument('-N', type=int, default=2, action="store", 
	help="Minimum number of reads supporting a UMI for it to be counted")
#parser.add_argument('-p', type=int, default=2, action="store",
#	help="Number of processors")

resample = parser.add_argument_group("resampling options")
resample.add_argument("-s", "--subsample", default=0, type=float, 
	help="Subsample reads to ratio of total reads (float) or to set value (int)")
resample.add_argument("-b", "--bootstrap", action="store_const", const=True, default=False, 
	help="Bootstrap quantification at full read depth")
resample.add_argument('--iter', type=int, default=42, 
	help="Number of iterations for resampling")

seal = parser.add_argument_group("Seal options")
seal.add_argument("-r", "--ref", type=str, default=None, 
	help="Reference to map FASTQ reads against")
seal.add_argument("-k", type=int, default=31, action="store",
	help="Kmer length for Seal (if you're processing reads and not seal output)")
seal.add_argument("-d", type=int, default=0, action="store",
	help="Set the 'hdist' parameter for Seal (# mismatches)")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

