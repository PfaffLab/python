#!/usr/bin/python
#==============================================================================
# scrna-prep-10x-reads.py
#
# Shawn Driscoll
# 20180124
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Translate 10x paired-end reads into a single-end file with the cell and umi
# barcodes incorporated into the read name. 
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser, basename, dirname
from collections import defaultdict
from time import localtime
from Basics import messages as ms
import gzip
from os import system
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock


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


	# variables
	left = []
	right = []

	# figure out what's up

	if args.a is not None or args.b is not None:
		if not (args.a is not None and args.b is not None):
			ms.error_message("You must specify both -a and -b if running a single sample.")
			return 1

	if args.a is not None and args.b is not None:
		# single sample
		left = [args.a]
		right = [args.b]

	elif args.f is not None:
		if not isfile(args.f):
			ms.error_message("Batch file does not exist")
			return 1

		# we have a batch file
		with open(args.f, "r") as fin:

			for szl in fin:
				aln = szl.strip().split("\t")
				left.append(aln[0])
				right.append(aln[1])

	if len(left) < 1 or len(right) < 1:
		ms.error_message("No samples to process!")
		return 1

	rres = core(left, right, args)
	
	return rres

def core(left, right, args):

	##
	## create queue and child process for compressing the fastq files
	##
	
	bc_len = args.barcode_length
	umi_len = args.umi_length
	
	tasks = JoinableQueue()
	p = Process(target=gz_worker, args=(tasks,))
	p.daemon = True
	p.start()

	##
	## input files are paired from the sequencer so we just have to read through them
	## and write them back out

	for i in range(len(left)):

		m1 = left[i]
		m2 = right[i]

		# pick apart the name of the second file to build the output name
		base = basename(m2)
		path = dirname(m2)

		base_parts = base.split(".")
		stub = base_parts[0]

		outfile = "{}_prepped.fastq".format(stub)

		if outfile==m1 or outfile==m2:
			ms.error_message("Output file path matches input file path. WTF? {}".format(outfile))
			sys.exit(1)

		if isfile(outfile):
			ms.warning_message("output file exists. overwriting. {}".format(outfile))

		fout = open(outfile, "w")

		try:
			ms.message("Processing {}".format(m1))
			with open_reads(m1) as fin1, open_reads(m2) as fin2:

				nidx = 0
				nreads = 0
				lread = []
				for szl2 in fin2:

					if (nreads % 1000000) == 0:
						ms.progress_message("Parsed {} reads".format(nreads))

					nidx += 1
					if nidx == 1:
						# read name line
						rname = szl2.strip().split()
						rname = rname[0]
						# read two lines from the barcode file
						szl1 = fin1.readline()
						szl1 = fin1.readline().strip()
						# this is the barcode so we can pick it apart. I'm going to put the cell barcode 
						# at the front of the read so that I can maybe leverage samtools sort to sort
						# barcodes together for me prior to parsing cells out
						rname_tmp = re.sub("^\@", "", rname)
						rname = "@{}:{}".format(szl1[0:bc_len], rname_tmp)
						# 20180226
						# moved the cell barcode to the front of the read name so we only
						# need to write the umi at the end and not both
						#rname += ":{}:{}".format(szl1[0:16], szl1[16:len(szl1)])
						rname += ":{}".format(szl1[bc_len:(bc_len+umi_len)])
						lread.append(rname + "\n")
						# read the remaining lines for this read from the barcodes file
						szl1 = fin1.readline()
						szl1 = fin1.readline()
					elif nidx < 4:
						lread.append(szl2)

					if nidx==4:
						# finished with read
						lread.append(szl2)
						fout.write("".join(lread))
						nidx = 0
						lread = []
						nreads += 1

			ms.progress_message("Parsed {} reads".format(nreads), last=True)

		except:
			fout.close()
			sys.exit(1)

		fout.close()

		#ms.message("compressing {}".format(outfile))
		#system("gzip -f {}".format(outfile))
		tasks.put(outfile)
	
	tasks.put(None)
	ms.message("Waiting for gzip compression to complete.")
	tasks.join()
	p.join()
	
	# done

	return 0

##
# child process function to gzip compress the fastq files
def gz_worker(qin):
	
	while True:
		item = qin.get()
		if item is None:
			qin.task_done()
			break
		
		
		if isfile(item):
			cmd = "gzip -f {}".format(item)
			system(cmd)
		
		qin.task_done()
	
	return
		


def open_reads(fname):

	if re.search("\.gz$", fname):
		return gzip.open(fname, "r")

	return open(fname, "r")

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Combine the 10x cell and umi barcodes into the read names of the second mates and output the second mate FASTQ.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('-a', type=str, default=None, help="First mate file (the barcodes file)")
parser.add_argument('-b', type=str, default=None, help="Second mate file (the actual reads)")
parser.add_argument('-f', type=str, default=None, 
	help="A two column text file with first and second mate files, first and second column respectivly, to process in batch.")

parser.add_argument('--barcode-length', type=int, default=16, help="Legth of cell barcode")
parser.add_argument('--umi-length', type=int, default=10, help="Length of UMI")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

