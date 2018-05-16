#!/usr/bin/python
#==============================================================================
# scrna-filter.py
#
# Shawn Driscoll
# 20171206
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Filters SCRNA reads. We can filter out reads that do not match anything in 
# a fasta reference OR filter out reads that match the reference. We can also 
# filter by UMI to exclude reads associated with UMIs that have frequency 
# lower than a specified amount.
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
import hashlib

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

	err = False

	##
	## confirm all files exist. check for target file.
	if not isfile(args.fin):
		error_message("Input file {} does not exist".format(args.fin))
		err = True

	if args.filter_match is not None:
		if not isfile(args.filter_match):
			error_message("Reference file {} does not exist".format(args.filter_match))
			err = True

	if args.filter_nonmatch is not None:
		if not isfile(args.filter_nonmatch):
			error_message("Reference file {} does not exist".format(args.filter_nonmatch))
			err = True

	if args.filter_nonmatch is not None and args.filter_match is not None:
		error_message("Please only specify one of --filter-match or --filter-nonmatch")
		err = True

	if err:
		return 1

	##
	## all good?
	##

	filtered = False

	ftype = get_filetype(args.fin)

	if not (ftype[0] == "fastq" or ftype[0] == "fasta"):
		error_message("Unable to detect file type for input reads")
		return 1

	tmpname = hashlib.md5(args.fin).hexdigest()
	if ftype[0]=="fastq":
		tmpname += ".fastq"
	else:
		tmpname += ".fasta"

	# setup output file name
	if ftype[1]:
		# file is gzipped
		tmp = (args.fin).split(".")
		n = len(tmp)
		tmpHat = tmp[0:(n-2)]
		tmpHat.append("filtered")
		tmpHat += tmp[(n-2):(n-1)]
		outfile = ".".join(tmpHat)
	else:
		# file is not gzipped
		tmp = (args.fin).split(".")
		n = len(tmp)
		tmpHat = tmp[0:(n-1)]
		tmpHat.append("filtered")
		tmpHat += tmp[(n-1):n]
		outfile = ".".join(tmpHat)

	
	if args.filter_nonmatch is not None:
		message("filtering out reads that match reference")
		rres = bbduk_filter(args.fin, outfile, args.filter_nonmatch, False, args.k)
		filtered = True

	if args.filter_match is not None:
		message("filtering out reads that do not match reference")
		rres = bbduk_filter(args.fin, outfile, args.filter_match, True, args.k)
		filtered = True

	if args.min_umi_freq > 1:
		message("filtering reads by UMI count")

		if filtered:
			# previous filtering was run so input is different
			cmd = "mv {} {}".format(outfile, tmpname)
			system(cmd)
			rres = umi_filter(tmpname, outfile, args.min_umi_freq, ftype[0], False)
			print rres
		else:
			rres = umi_filter(infile, outfile, args.min_umi_freq, ftype[0], ftype[1])
			print rres


	if ftype[1]:
		message("input was gzip'd so we're gzipping the output")
		cmd = "gzip {}".format(outfile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)

	if isfile(tmpname):
		message("removing temporary file")
		cmd = "rm {}".format(tmpname)
		system(cmd)


	return 0

def get_filetype(fname):

	tmp = str.lower(fname).split(".")
	n = len(tmp)

	gz = tmp[-1] == "gz"

	if gz:
		ext = tmp[-2]
	else:
		ext = tmp[-1]

	if ext == "fasta" or ext == "fa":
		fmt = "fasta"
	elif ext == "fastq" or ext == "fq":
		fmt = "fastq"

	return [fmt, gz]

##
# this function returns the cell barcode and umi barcode from a read name string
def parse_barcodes(sz):

	tmp = sz.split(":")
	umi = tmp[-1]
	cell = tmp[-2]

	return [cell, umi]

##
# runs bbduk to filter reads against a reference. select is a boolean. if true then
# the function keeps matched reads otherwise it keeps unmatched reads.
def bbduk_filter(infile, outfile, ref, select, k):

	cmd = "bbduk.sh in={} ref={} k={} hdist=0 mm=f overwrite=t".format(infile, ref, k)
	if select:
		cmd += " outm={}".format(outfile)
	else:
		cmd += " out={}".format(outfile)

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

parser.add_argument('fin', type=str, help="Input fasta/fastq gz or not")

parser.add_argument("--filter-match", type=str, default=None, 
	help="FASTA reference to filter reads against. Matched reads are retained.")

parser.add_argument("--filter-nonmatch", type=str, default=None, 
	help="FASTA reference to filter reads against. Matched reads are discarded.")

parser.add_argument("-k", type=int, default=31, 
	help="Kmer size for read filtering against a reference.")

parser.add_argument('--min-umi-freq', type=int, default=1, 
	help="Minimum frequency of a UMI within a cell for output to file.")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

