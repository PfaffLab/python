#!/usr/bin/python
#==============================================================================
# demux-cellranger-aligned.py
#
# Shawn Driscoll
# 20170914
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
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime, time
import subprocess as sp
from os import system
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock

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
	bam_file = False

	bc = {}
	bc_umi = {}
	num_bc = 0
	offset = 0
	lnum = 0
	sz_umi = ""
	umi_file= ""
	# string to capture file summary table that's used with kallisto pseudo -b
	sz_table = "#id\tumiFile\tcellFile\n" 
	dumi = None

	filter_umi = args.min_umi_freq > 1

	if filter_umi:
		sys.stderr.write("NOTE: UMI's with frequency less than {} will be dropped\n".format(args.min_umi_freq))

	file_queue = JoinableQueue()
	child_p = None

	t0 = time()
	#
	# open input reads file. may be either sam or bam. reads are filtered while being 
	# read by samtools to remove secondary alignments. STAR is good about this so you 
	# do end up with distinct reads only
	if re.search("\.sam$", args.fin):
		p0 = sp.Popen("samtools view -hSF 0x100 {}".format(args.fin).split(), stdout=sp.PIPE)
	elif re.search("\.bam$", args.fin):
		bam_file = True
		p0 = sp.Popen("samtools view -hF 0x100 {}".format(args.fin).split(), stdout=sp.PIPE)
	else:
		message("Unknown input file format")
		return(1)
	
	t0 = time()
	#
	# parse the alignments. in this loop we only extract the cell barcode and the umi
	# plus record the file position offsets for barcodes. the dict that is built
	# is indexed by the barcodes and each element contains a list of file offsets for 
	# reads that came from that barcode. we also get all of the distinct umis collected 
	# per barcode in this loop in order to estimate the actual cell count before
	# writing all of the read files out to disk
	message('Parsing alignments from cellranger')
	with p0.stdout as fin, open("temp.sam", "w") as fout:

		for szl in fin:
			# write sam files that pass filter out to disk. this is the file that 
			# will be used to access the reads later
			fout.write(szl)
			if szl[0]=="@":
				# write header lines out to the file as well. record offset.
				offset += len(szl)
				continue

			# count lines and produce progress message so we know this thing is 
			# running
			lnum += 1
			if lnum % 1000000 == 0:
				progress_message("read {} lines".format(lnum))

			# get the cell id. this is the cell barcode with '-1' appended to it for a normal 
			# single sample experiment. if multiple experiments were aggregated these tags will 
			# be distinct per sample. this way you get output files for each original sample. 
			line_bc = parse_barcode(szl)
			if line_bc is None:
				# sometimes there is no cell id. not sure how that happens but whatever. skip
				# this read.
				offset += len(szl)
				continue

			if line_bc not in bc:
				# increment barcode count
				num_bc += 1
				# init a list for this barcode
				bc[line_bc] = []
				# init a dict for the barcode to track umis
				bc_umi[line_bc] = defaultdict(int)

			# append line offset to this barcode's list
			bc[line_bc].append(offset)
			# get the umi and add it to this barcode's dict
			umi = parse_umi(szl)
			bc_umi[line_bc][umi] += 1

			# update offset to the next line
			offset += len(szl)

	# final progress message and total time of parsing
	progress_message("read {} lines".format(lnum), last=True)
	sys.stderr.write("{} sec\n".format(time()-t0))

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
	with open("barcode_umi_counts.txt", "w") as fout:
		bc_umi_counts = []
		for lbc in bc.keys():
			num_umi = len(bc_umi[lbc].keys())
			bc_umi_counts.append([lbc, num_umi])
			# write the cell id, distinct umi count and total read count to file
			fout.write("\t".join(map(str, [lbc, num_umi, len(bc[lbc])])))
			fout.write("\n")

	# sort by umi count in descending order and threshold
	bc_umi_counts.sort(key=lambda x: x[1], reverse=True)
	exp_cells = int(math.floor(args.exp_cells*0.01 - 1))
	i = 0
	while True:
		if bc_umi_counts[i][1] < bc_umi_counts[exp_cells][1]*1.0/10:
			break
		i += 1
	# number of actual cells is 'i' because 'i' is incremented before 
	# checking if the umi count passes the threshold. i-1 is the index
	# of the last cell we would accept
	num_cells = i

	#
	# let user know what's up
	sys.stderr.write("{} sec\n".format(time()-t0))
	sys.stderr.write("Total distinct barcodes:  {}\n".format(num_bc))
	sys.stderr.write("Cell number estimate:     {}\n".format(num_cells))

	if args.force_cells is not None:
		# change number of cells to either the total barcodes or the 
		# value provided by the user, whichever is smaller
		num_cells = min([args.force_cells, num_bc])
		sys.stderr.write("Forced cell output:       {}\n".format(num_cells))

	t0 = time()
	message("Writing cell files")
	fin = open("temp.sam", "r")

	# start child process for compressing reads files. we'll use a 
	# queue for that so it can just deal with files as they are finished
	# being written and run independently of the main process.
	message("Launching child process for compressing read files")
	child_p = Process(target=compress_reads, args=(file_queue,))
	child_p.daemon = True
	child_p.start()

	# write individual cell files
	i = 0
	sz_umi = ""
	while i < num_cells:
		# get barcode
		lbc = bc_umi_counts[i][0]
		# start output strings
		szout = ""
		sz_umi = ""
		# setup output file name
		if args.fasta:
			cell_file = "{}.fasta".format(lbc)
		elif args.fastq:
			cell_file = "{}.fastq".format(lbc)
		else:
			cell_file = "{}.sam".format(lbc)
		
		umi_file = "{}.umi".format(lbc)

		# update user on progress
		progress_message("Writing {} - {}/{} ({} reads)".format(cell_file, i+1, num_cells, len(bc[lbc])))

		if filter_umi:
			##
			## build a dict with umi's and their frequency within this barcode.
			##
			dumi = defaultdict(int)
			for offset in bc[lbc]:
				fin.seek(offset)
				szl = fin.readline().strip()
				umi = parse_umi(szl)
				dumi[umi] += 1

		# loop through line offsets for this barcode
		for offset in bc[lbc]:
			fin.seek(offset)
			szl = fin.readline().strip()
			umi = parse_umi(szl)
			if filter_umi:
				if dumi[umi] < args.min_umi_freq:
					# drop this one
					continue

			if umi is not None:
				
				if not args.no_umi_file:
					sz_umi += "{}\n".format(umi)
				
				aln = szl.split("\t")
				# append the cell barcode and the umi barcode to the read name
				aln[0] += ":{}:{}".format(strip_cellid(lbc), umi)

				if args.fasta:
					r = aln[9]
					if (int(aln[1]) & 0x10) != 0:
						temp = revcomp(r)
						r = temp
					szl = ">{}\n{}\n".format(aln[0], r)
				elif args.fastq:
					r = aln[9]
					qual = aln[10]
					if (int(aln[1]) & 0x10) != 0:
						r = revcomp(r)
						temp = list(qual)
						temp.reverse()
						qual = "".join(temp)
					szl = "@{}\n{}\n+\n{}\n".format(aln[0], r, qual)
				else:
					szl = "\t".join(aln) + "\n"

				szout += szl
		
		with open(cell_file, "w") as fout:
			fout.write(szout)
		
		if not args.no_umi_file:
			with open(umi_file, "w") as fout:
				fout.write(sz_umi)

		file_queue.put(cell_file)

		sz_table += "{}\t{}\t{}.gz\n".format(lbc, umi_file, cell_file)

		i += 1

	fin.close()

	sys.stderr.write("\n")
	sys.stderr.write("{} sec\n".format(time()-t0))

	with open("sample_file_table.tsv", "w") as fout:
		fout.write(sz_table)
	
	sys.stderr.write("Waiting for child process to finish compressing files\n")
	file_queue.put(None)
	file_queue.join()
	child_p.join()

	system("rm -f temp.sam")
	message("finished!")

	return 0

##
## this code extracts the cell id just in case we are dealing with an aggregated 
## set of samples. this function drops off the '-[0-9]' suffix so we just have the 
## barcode.
def strip_cellid(cellid):
	tmp = cellid.split("-")
	return(tmp[1])

def compress_reads(tasks):

	while True:
		item = tasks.get()
		if item is None:
			tasks.task_done()
			break


		cmd = "gzip -f {}".format(item)
		system(cmd)
		tasks.task_done()


def parse_barcode(sz):
	r = re.search("CB:Z:([^\s]+)", sz)
	if r:
		return r.group(1)

	return None

def parse_umi(sz):
	r = re.search("UB:Z:([^\s]+)", sz)
	if r:
		return r.group(1)

	return None

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

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('fin', type=str, help="Input file (SAM or BAM)")

parser.add_argument('--no-umi-file', default=False, action="store_const", const=True, 
	help="Do not write separate files with the UMI barcodes. These files are used if you use kallisto 'pseduo' to quantify.")

parser.add_argument('--exp-cells', type=int, default=3000, 
	help="Expected number of cells.")

parser.add_argument("--force-cells", type=int, default=None, 
	help="Force this many cells to be extracted. Sorted from top umi count.")

parser.add_argument('--min-umi-freq', type=int, default=2, 
	help="Minimum frequency of a UMI within a cell for output to file.")

mgroup = parser.add_mutually_exclusive_group(required=False)

mgroup.add_argument("--fasta", action="store_const", const=True, default=False, 
	help="Write reads in FASTA format")
mgroup.add_argument("--fastq", action="store_const", const=True, default=False, 
	help="Write reads in FASTA format")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

