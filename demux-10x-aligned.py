#!/usr/bin/python
#==============================================================================
# demux-cellranger-aligned.py
#
# Shawn Driscoll
# 20180207
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses the "bam" alignment output file from cellranger to produce individual
# cell read files for manual analysis. Output may be either in bam, fasta or 
# fastq format.  Writing FASTA is significantly faster than the other formats
# and if your downstream analysis does not require the invidiual base quals
# (aligners like STAR ignore these as well as quantification programs like kallisto
# and sailfish/salmon) then it may be desireable to export in this format. 
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from collections import defaultdict
from time import localtime, time
import subprocess as sp
import os, os.path
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock
from Basics import messages as ms
import cPickle as pickle
import random

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
BC_PICKLE = "bc_index.pkl"
BC_READCOUNT = "bc_readcount.pkl"
BC_UMI_PICKLE = "bc_umi_index.pkl"

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	bam_file = False

	bc = {}
	bc_readcount = defaultdict(int)
	bc_umi = {}
	num_bc = 0
	offset = 0
	lnum = 0
	sz_umi = ""
	umi_file= ""
	# string to capture file summary table that's used with kallisto pseudo -b
	sz_table = "#id\tumiFile\tcellFile\n"
	dumi = None
	bam_flag = False
	sam_header = ""
	quant_mode = False

	##
	## check for output folder
	##
	if not os.path.isdir(args.outpath):
		ms.message("Creating output folder {}".format(args.outpath))
		os.mkdir(args.outpath)

	file_queue = JoinableQueue()
	p = None
	pool = []

	if args.R is not None:
		if not os.path.isfile(args.R):
			ms.error_message("Supplied annotation file does not exist ({})".format(args.R))
			return 1
		else:
			quant_mode = True

	##
	## figure out if we have a bam as input. if so we have to convert it to sam for indexing
	##
	if re.search("\.bam$", args.fin):
		# send the sam file into the output folder
		sam_name = args.outpath + "/" + os.path.basename(re.sub("\.bam$", ".sam", args.fin))
		
		if not os.path.isfile(sam_name):
			# need to convert alignments to sam
			bam_flag = True
			cmd = "samtools view -h {} > {}".format(args.fin, sam_name)
			t0 = time()
			message("Temporarily converting BAM to SAM format")
			rres = runcmd(cmd)
			if rres[0] != 0:
				sys.stderr.write("Error: samtools exited with non-zero exit status!\n")
				return 1
				
			sys.stderr.write("{} sec\n".format(time()-t0))
				
	else:
		sam_name = args.fin 
	
	
	##
	## we need to index all barcodes and track umi per barcode. if these pickle 
	## files exist we can use them
	##
	
	bc_pkl = args.outpath + "/" + BC_PICKLE
	bc_umi_pkl = args.outpath + "/" + BC_UMI_PICKLE 
	bc_readcount_pkl = args.outpath + "/" + BC_READCOUNT
	sam_header_pkl = args.outpath + "/sam_header.pkl"
	
	if os.path.isfile(bc_pkl) and os.path.isfile(bc_umi_pkl) and os.path.isfile(bc_readcount_pkl):
		##
		# load indexes from pickles
		ms.message("Loading existing barcode and umi indexes from output folder")
		t0 = time()
		bc = pickle.load(open(bc_pkl, "rb"))
		bc_umi = pickle.load(open(bc_umi_pkl, "rb"))
		bc_readcount = pickle.load(open(bc_readcount_pkl, "rb"))
		sam_header = pickle.load(open(sam_header_pkl, "rb"))
		num_bc = len(bc.keys())
		ms.time_diff(t0)
	
	else:
		# we have to index 
	
		#
		# parse the alignments. in this loop we only extract the cell barcode and the umi
		# plus record the file position offsets for barcodes. the dict that is built
		# is indexed by the barcodes and each element contains a list of file offsets for 
		# reads that came from that barcode. we also get all of the distinct umis collected 
		# per barcode in this loop in order to estimate the actual cell count before
		# writing all of the read files out to disk
		message('Indexing cell barcodes from alignments and counting raw UMI.')
		t0 = time()
		with open(sam_name, "r") as fin:
	
			for szl in fin:
				if szl[0]=="@":
					# append header line to header string
					sam_header += szl
					offset += len(szl)
					continue
	
				# count lines and produce progress message so we know this thing is 
				# running
				lnum += 1
				if lnum % 1000000 == 0:
					progress_message("read {} lines".format(lnum))
	
				# fetch the cell barcode from the read name
				line_bc = parse_barcode(szl)
	
				if line_bc not in bc:
					# first encounter with this barcode
					num_bc += 1
					# init a list for this barcode's line offsets within this sam file
					bc[line_bc] = []
					# init a dict for the barcode to track umis
					bc_umi[line_bc] = defaultdict(int)
	
				# append line offset to this barcode's list
				bc[line_bc].append(offset)
				# get the umi and add it to this barcode's dict IF this is not a
				# secondary alignment
				aln = szl.split("\t")
				
				if (int(aln[1]) & 0x100) == 0:
					# not a secondary alignment. track it.
					umi = parse_umi(szl)
					bc_umi[line_bc][umi] += 1
				
				if ((int(aln[1]) & 0x4) == 0) and((int(aln[1]) & 0x100) == 0):
					# this read is aligned and is a primary alignment so we can count this one
					# into this barcode's aligned read count
					bc_readcount[line_bc] += 1
	
				# update offset to the next line
				offset += len(szl)
	
		# final progress message and total time of parsing
		progress_message("read {} lines".format(lnum), last=True)
		sys.stderr.write("{} sec\n".format(time()-t0))
	
		t0 = time()	
		
		if not args.no_pickles:
			ms.message("saving indexes to disk")
			pickle.dump(bc, open(bc_pkl, "wb"))
			pickle.dump(bc_umi, open(bc_umi_pkl, "wb"))
			pickle.dump(bc_readcount, open(bc_readcount_pkl, "wb"))
			pickle.dump(sam_header, open(sam_header_pkl, "wb"))
			ms.time_diff(t0)

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
	message("Determining cell count")
	
	#
	# write a file that will contain the cell id, umi count and read count
	# for each cell id. might be informative...who knows.
	with open("{}/barcode_umi_counts.txt".format(args.outpath), "w") as fout:
		bc_umi_counts = []
		
		fout.write("barcode\tumi_count\tdistinct_reads\talignments\n")
		
		for lbc in bc.keys():
			num_umi = len(bc_umi[lbc].keys())
			bc_umi_counts.append([lbc, num_umi])
			# write the cell id, distinct umi count and total read count to file
			fout.write("\t".join(map(str, [lbc, num_umi, bc_readcount[lbc], len(bc[lbc])])))
			fout.write("\n")
	
	#
	# sort by umi count in descending order and threshold
	bc_umi_counts.sort(key=lambda x: x[1], reverse=True)
	exp_cells = int(math.floor(args.exp_cells*0.01 - 1))
	
	num_reads = 0
	num_umi = 0
	
	i = 0
	while True:
		if bc_umi_counts[i][1] < bc_umi_counts[exp_cells][1]*1.0/10:
			break
		
		# count umi and count distinct reads
		lbc = bc_umi_counts[i][0]
		num_reads += bc_readcount[lbc]
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
		fout.write("median_reads_per_cell\t{}\n".format(bc_readcount[median_lbc]))

	#
	# let user know what's up
	sys.stderr.write("{} sec\n".format(time()-t0))
	sys.stderr.write("Total distinct barcodes:  {}\n".format(num_bc))
	sys.stderr.write("Cell number estimate:     {}\n".format(num_cells))

	if args.estimate_only:
		if bam_flag:
			# input was BAM so we can dump the converted file. just putting in
			# some logic to be certain that the original file is not deleted.
			if os.path.isfile(args.fin) and os.path.isfile(sam_name) and (sam_name != args.fin):
				os.unlink(sam_name)

		ms.message("Done.")
		return 0

	if args.force_cells is not None:
		# change number of cells to either the total barcodes or the 
		# value provided by the user, whichever is smaller
		num_cells = min([args.force_cells, num_bc])
		sys.stderr.write("Forced cell output:       {}\n".format(num_cells))

	
	t0 = time()
	
	message("Parsing individual detected cell alignments out to individual files")

			
	if quant_mode:
		# start quantificaion child processes for parsed sam files
		for i in range(args.p):
			p = Process(target=quantification_worker, args=(file_queue, args, ))
			p.daemon = True
			p.start()
			pool.append(p)

	else:
		# start child process for sam to bam conversion 
		for i in range(args.p):
			p = Process(target=compress_reads, args=(file_queue,))
			p.daemon = True
			p.start()
			pool.append(p)

	fin = open(sam_name, "r")

	# write individual cell files
	i = 0
	sz_umi = ""
	while i < num_cells:
		# get barcode
		lbc = bc_umi_counts[i][0]
		# start output strings
		szout = sam_header
		#sz_umi = ""
		# setup output file name
		cell_file = "{}/{}.sam".format(args.outpath, lbc)
		#umi_file = "{}.umi".format(lbc)

		# update user on progress
		progress_message("Writing {} - {}/{} ({} reads)".format(cell_file, i+1, num_cells, len(bc[lbc])))

		if args.samplerate < 1 and args.samplerate > 0:
			
			##
			# to subsample we have to run through all read offsets for this cell and index the reads
			# then take a subset of them to write out to disk. I have to do this because the 
			# alignment file contains secondary alignments which have to be collapsed by 
			# read name prior to the subsampling.
			read_index = defaultdict(list)
			for offset in bc[lbc]:
				fin.seek(offset)
				aln = fin.readline().strip().split("\t")
				rname = aln[0]
				read_index[rname].append(offset)
			
			#
			# now by looping through distinct reads we can dump out only those that are at the specified rate
			for rname in read_index.keys():
				if random.random() > args.samplerate:
					continue
				
				# dump this read
				for offset in read_index[rname]:
					fin.seek(offset)
					szout += fin.readline()
						
		else:

			# loop through line offsets for this barcode and append lines to the output string
			for offset in bc[lbc]:
				fin.seek(offset)
				szout += fin.readline()
		
		# write the file
		with open(cell_file, "w") as fout:
			fout.write(szout)
		
		# send the file off for bam compression
		file_queue.put(cell_file)

		i += 1

	fin.close()

	sys.stderr.write("\n")
	sys.stderr.write("{} sec\n".format(time()-t0))
	
	if bam_flag:
		# input was BAM so we can dump the converted file. just putting in
		# some logic to be certain that the original file is not deleted.
		if os.path.isfile(args.fin) and os.path.isfile(sam_name) and (sam_name != args.fin):
			os.unlink(sam_name)
			
	sys.stderr.write("Waiting for child process to finish compressing files\n")
	
	for p in pool:
		file_queue.put(None)
	file_queue.join()
	
	for p in pool:
		p.join()

	message("finished!")

	return 0

def runcmd(cmd, verbose=True):
	
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	
	rres = os.system(cmd)
	
	signal_num = rres & 127
	exit_status = rres >> 8
	
	return (exit_status, signal_num)

def compress_reads(tasks):

	while True:
		item = tasks.get()
		if item is None:
			tasks.task_done()
			break

		bam_name = re.sub("\.sam$", ".bam", item)
		if os.path.isfile(item):
			cmd = "samtools view -bS -o {} {}".format(bam_name, item)
			rres = runcmd(cmd, verbose=False)
			if rres[0] == 0 and os.path.isfile(bam_name):
				os.unlink(item)
		
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
		

def parse_barcode(sz):
	# break up the line and then the read name to extract the cell barcode
	# since the read names are "namepartsandjunk:cell_barcode:umi"
	aln = sz.split()
	tmp = aln[0].split(":")
	n = len(tmp)
	return tmp[n-2]	

def parse_umi(sz):
	# break up the line and then the read name to extract the umi
	# since the read names are "namepartsandjunk:cell_barcode:umi"
	aln = sz.split()
	tmp = aln[0].split(":")
	n = len(tmp)
	return tmp[n-1]	

def revcomp(sz):
	temp = list(sz.upper())
	change = { "A":"T", "G":"C", "T":"A", "C":"G", "N":"N" }
	flip = [change[k] for k in temp]
	flip.reverse()
	return "".join(flip)


def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0


def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

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

parser = argparse.ArgumentParser(description="Assuming 10x reads mapped to a genome in a BAM file, this program can estimate number of detected cells and extract their alignments out to individual BAM files OR quantify UMI expression per cell.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('outpath', type=str, help="Folder for output. It will be created if it does not exist already.")
parser.add_argument('fin', type=str, help="Input file (SAM or BAM)")

parser.add_argument('--estimate-only', action="store_const", const=True, default=False, 
	help="Run up to estimation of detected cells and dump information to output folder.")

parser.add_argument('-p', type=int, default=1, help="Number of child processes for BAM compression or quantification.")

parser.add_argument('--no-pickles', action="store_const", const=True, default=False, 
	help="Do not save the barcode and UMI index dicts as pickles")

parser.add_argument('--exp-cells', type=int, default=3000, 
	help="Expected number of cells.")

parser.add_argument("--force-cells", type=int, default=None, 
	help="Force this many cells to be extracted. Sorted from top umi count.")

parser.add_argument("--samplerate", type=float, default=1, 
	help="Subsample distinct reads written out to cell files at this rate. Useful for when you need to normalize reads per cell between separate libraries.")

parser.add_argument('-R', type=str, default=None, help="RefFlat annotation for UMI quantification.")
parser.add_argument('-q', type=int, default=0, help="Minimum MAPQ for running quantification step.")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

