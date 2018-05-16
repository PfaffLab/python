#!/usr/bin/python
#==============================================================================
# demux-10x-sorted.py
#
# Shawn Driscoll
# 20180226
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses a name sorted BAM to estimate cell count and identify cells. Input 
# should be a name-sorted BAM file of unaligned reads with the cell barcodes
# at the front of the read names (use scrna-prep-10x-reads.py), convert to 
# BAM and then name sort with 'samtools sort -n ...' 
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from collections import defaultdict
import pysam as ps
from time import localtime, time
import subprocess as sp
import os, os.path
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock
from Basics import messages as ms
#import cPickle as pickle
import random
import gzip

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

HOME = os.path.expanduser("~")
PYTHON_PATH = "{}/coding/python".format(HOME)

#==============================================================================
# main
#==============================================================================

def main(args):
	
	# check input file
	if not os.path.isfile(args.bam):
		ms.error_message("Input file does not exist")
		return 1
	
	if not re.search("\.bam$", args.bam):
		ms.error_message("Input file should be a BAM file.")
		return 1
		
	
	##
	## check for output folder
	##
	if not os.path.isdir(args.outpath):
		ms.message("Creating output folder {}".format(args.outpath))
		os.mkdir(args.outpath)

	rres = core(args)
	
	return rres
	
	
	
def core(args):

	# variables
	
	# barcode dict to count reads per barcode
	bc = defaultdict(int)
	# barcode dict to track distinct umi and count them
	bc_umi = {}
	bc_keep = set()
	num_bc = 0
	lnum = 0
	lbc = ""
	lbc_last = ""
	
	# for output conversion
	file_queue = JoinableQueue()
	p = None
	pool = []
		
	##
	## we need to index all barcodes and track umi per barcode. if these pickle 
	## files exist we can use them
	##
	
	# we have to index 

	#
	# parse the alignments. in this loop we only extract the cell barcode and the umi
	# plus record the file position offsets for barcodes. the dict that is built
	# is indexed by the barcodes and each element contains a list of file offsets for 
	# reads that came from that barcode. we also get all of the distinct umis collected 
	# per barcode in this loop in order to estimate the actual cell count before
	# writing all of the read files out to disk
	ms.message('Counting per-barcode reads and UMI')
	t0 = time()
	with ps.AlignmentFile(args.bam, "rb", check_header=False, check_sq=False) as fin:

		for aln in fin:
			lnum += 1
			
			nparts = (aln.query_name).split(":")
			# barcode is first
			lbc = nparts[0]
			# umi is last
			umi = nparts[-1]
			
			bc[lbc] += 1
			
			if lbc not in bc_umi:
				bc_umi[lbc] = {}
			
			if umi not in bc_umi[lbc]:
				bc_umi[lbc][umi] = 0
			
			bc_umi[lbc][umi] += 1
			
			if (lnum % 1000000) == 0:
				ms.progress_message("parsed {} reads".format(lnum))
				

	# final progress message and total time of parsing
	ms.progress_message("parsed {} reads".format(lnum), last=True)
	sys.stderr.write("{} sec\n".format(time()-t0))

	t0 = time()	

	#
	# implement cell number detection per 10x.
	# here's what happens. you take the 'exp-cells' value (expected cells)
	# and multiply that by 0.01 to get an index. sort the barcodes and the 
	# barcode umi counts in descending order and jump to the index you just
	# calculated and then take that index's umi count. scale that count 
	# by 0.1. now you take as many cells, starting from the top of the umi
	# count sorted list, that have at least that many UMI.  that's literally 
	# how they do it. 
	#

	t0 = time()
	ms.message("Determining cell count")
	
	#
	# write a file that will contain the cell id, umi count and read count
	# for each cell id. might be informative...who knows.
	with open("{}/barcode_umi_counts.txt".format(args.outpath), "w") as fout:
		bc_umi_counts = []
		
		fout.write("barcode\tumi_count\tdistinct_reads\n")
		
		for lbc in bc.keys():
			num_umi = len(bc_umi[lbc].keys())
			bc_umi_counts.append([lbc, num_umi])
			# write the cell id, distinct umi count and total read count to file
			fout.write("\t".join(map(str, [lbc, num_umi, bc[lbc]])))
			fout.write("\n")
	
	#
	# sort by umi count in descending order and threshold
	bc_umi_counts.sort(key=lambda x: x[1], reverse=True)
	exp_cells = int(math.floor(args.exp_cells*0.01 - 1))
	
	num_reads = 0
	num_umi = 0
	num_bc = len(bc.keys())
	
	i = 0
	while True:
		# check if the current barcode is below threshold..
		if bc_umi_counts[i][1] < bc_umi_counts[exp_cells][1]*1.0/10:
			break
		
		# count umi and count distinct reads
		lbc = bc_umi_counts[i][0]
		# keep track of the barcodes that we will retain
		num_reads += bc[lbc]
		num_umi += len(bc_umi[lbc].keys())
		
		i += 1
	
	#
	# number of actual cells is 'i' because 'i' is incremented before 
	# checking if the umi count passes the threshold. i-1 is the index
	# of the last cell we would accept
	num_cells = i
	
	#
	# now we can generate a summary for the detected cells
	with open("{}/cell_summary.tsv".format(args.outpath), "w") as fout:

		fout.write("estimated_cells\t{}\n".format(num_cells))
		fout.write("total_reads\t{}\n".format(num_reads))
		fout.write("total_umi\t{}\n".format(num_umi))
		fout.write("reads_per_cell\t{}\n".format(num_reads * 1.0/num_cells))
		fout.write("umi_per_cell\t{}\n".format(num_umi * 1.0/num_cells))

		# find the median barcode and corresponding read count
		if num_cells % 2 == 0:
			# even count
			median_idx = num_cells/2
		else:
			median_idx = num_cells/2 + 1
		
		median_lbc = bc_umi_counts[median_idx][0]
		fout.write("median_reads_per_cell\t{}\n".format(bc[median_lbc]))

	#
	# let user know what's up
	sys.stderr.write("{} sec\n".format(time()-t0))
	sys.stderr.write("Total distinct barcodes:  {}\n".format(num_bc))
	sys.stderr.write("Cell number estimate:     {}\n".format(num_cells))

	if args.estimate_only:
		ms.message("Done.")
		return 0

	if args.force_cells is not None:
		# change number of cells to either the total barcodes or the 
		# value provided by the user, whichever is smaller
		num_cells = min([args.force_cells, num_bc])
		sys.stderr.write("Forced cell output:       {}\n".format(num_cells))

	t0 = time()

	# make set of the barcodes that we will keep
	bc_keep = set()
	for i in range(num_cells):
		bc_keep.add(bc_umi_counts[i][0])

	#
	# now we can dig back into the sorted bam file to export all of the individual cell lines
	
	# launch processes for gzip compression
	for i in range(args.p):
		p = Process(target=compress_reads, args=(file_queue,))
		p.daemon = True
		p.start()
		pool.append(p)
	
	with ps.AlignmentFile(args.bam, "rb", check_header=False, check_sq=False) as fin:
		
		szout = ""
		bc_out = 0
		for aln in fin:
			lnum += 1
			
			nparts = (aln.query_name).split(":")
			# barcode is first
			lbc = nparts[0]
			# umi is last
			#umi = nparts[-1]
			
			if lbc != lbc_last:
				if lbc_last in bc_keep:
					
					# write buffered data to file...
					bc_out += 1
					fname = "{}/{}.fastq".format(args.outpath, lbc_last)
					ms.progress_message("writing {} ({} of {})".format(fname, bc_out, num_cells))
					with open(fname, "w") as fout:
						fout.write(szout)
					
					file_queue.put(fname)
					
				szout = ""
			
			# keep it?
			if random.random() < args.samplerate:
				# passes sampling limit
				if lbc in bc_keep:
					# convert line to fasta and append it to the output string
					szout += fastq_from_aln(aln)
				
			lbc_last = lbc
			
	if lbc_last in bc_keep:
		bc_out += 1
		# write buffered data to file...
		fname = "{}/{}.fastq".format(args.outpath, lbc_last)
		ms.progress_message("writing {} ({} of {})".format(fname, bc_out, num_cells), last=True)
		with open(fname, "w") as fout:
			fout.write(szout)
		
		# put fname in the queue to be compressed
		file_queue.put(fname)
	
	for p in pool:
		file_queue.put(None)
	
	file_queue.join()
	
	for p in pool:
		p.join()

	ms.message("finished!")

	return 0

##
# make a fastq version of a SAM alignment. returns as a single string with
# all four lines separated by \n with an additional \n at the end.
# this function is slightly specialized for this file's specific goal.
def fastq_from_aln(aln):
	
	# reformat the read name to be name:barcode:umi
	tmp = (aln.query_name).split(":")
	lbc = tmp[0]
	umi = tmp[-1]
	n = len(tmp)
	
	r0 = ":".join(tmp[1:(n-2)])
	
	rname = "@{}:{}:{}".format(r0, lbc, umi)
	
	seq = aln.query_sequence
	qual = aln.qual

	sz = "\n".join([rname, seq, "+", qual]) + "\n"
	return(sz)

def runcmd(cmd, verbose=True):
	
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	
	rres = os.system(cmd)
	
	signal_num = rres & 127
	exit_status = rres >> 8
	
	return (exit_status, signal_num)

##
# worker function for compressing the fastq files with gzip.
def compress_reads(tasks):

	while True:
		item = tasks.get()
		if item is None:
			tasks.task_done()
			break

		rres = runcmd("gzip {}".format(item), verbose=False)
		#if rres[0] != 0:
			
		tasks.task_done()
	
	return 0

def quantification_worker(tasks, args):
	
	py_conv = "{}/bam-galign-to-talign-fast.py".format(PYTHON_PATH)
	py_quant = "{}/scrna-samhits.py".format(PYTHON_PATH)
	
	while True:
		item = tasks.get()
		if item is None:
			tasks.task_done()
			break
		
		if not os.path.isfile(item):
			tasks.task_done()
			continue
		
		# 'item' is the filename of SAM alignments for a single cell
		tsam_out = re.sub("\.sam", "_talign.sam", item)
		thits_out = re.sub("\.sam", ".umi_count", item)
		
		cmd = py_conv + " --filter-mapq 0 --filter-overlap 0.2 --f-stranded -o {} {} {} 2>/dev/null".format(tsam_out, args.R, item)
		rres = runcmd(cmd, verbose=False)
		if rres[0] == 0:
			# normal exit. good.
			cmd = py_quant + " -q {} --quiet {} {} > {} 2>/dev/null".format(args.q, args.R, tsam_out, thits_out)
			rres = runcmd(cmd, verbose=False)
		
		if os.path.isfile(tsam_out):
			os.unlink(tsam_out)
		
		# 
		# remove the sam file?
		if os.path.isfile(item):
			os.unlink(item)
		
		tasks.task_done()
	
	return 0			
		
def revcomp(sz):
	temp = list(sz.upper())
	change = { "A":"T", "G":"C", "T":"A", "C":"G", "N":"N" }
	flip = [change[k] for k in temp]
	flip.reverse()
	return "".join(flip)


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Unaligned, read name sorted with cell barcode in FRONT of name, BAM as input. Output will be a couple text file reports and each cell's reads as a FASTQ.gz file.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('outpath', type=str, help="Folder for output. It will be created if it does not exist already.")
parser.add_argument('bam', type=str, help="Input file (SAM or BAM)")

parser.add_argument('--estimate-only', action="store_const", const=True, default=False, 
	help="Run up to estimation of detected cells and dump information to output folder.")

parser.add_argument('-p', type=int, default=1, help="Number of child processes for gzip compression.")

parser.add_argument('--exp-cells', type=int, default=3000, 
	help="Expected number of cells.")

parser.add_argument("--force-cells", type=int, default=None, 
	help="Force this many cells to be extracted. Sorted from top umi count.")

parser.add_argument("--samplerate", type=float, default=1, 
	help="Subsample distinct reads written out to cell files at this rate. Useful for when you need to normalize reads per cell between separate libraries.")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

