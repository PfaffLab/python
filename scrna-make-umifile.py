#!/usr/bin/python
#==============================================================================
# scrna-make-umifile.py
#
# Shawn Driscoll
# 20171207
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Extracts the UMI barcode from reads and writes them out to a plain text file
# with one line per read. Needed if you're going to run kallisto pseudo. The
# UMI barcode should be the last thing attached to the read name and must be 
# prefixed by a ':'.
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
	## confirm input file exists
	##
	
	if not isfile(args.fin):
		error_message("Input file {} does not exist".format(args.fin))
		err = True

	if err:
		return 1

	##
	## all good?
	##

	ftype = get_filetype(args.fin)

	if not (ftype[0] == "fastq" or ftype[0] == "fasta"):
		error_message("Unable to detect file type for input reads")
		return 1

	# setup output file name
	if ftype[1]:
		# file is gzipped
		tmp = (args.fin).split(".")
		n = len(tmp)
		tmpHat = tmp[0:(n-2)]
		tmpHat.append("umi")
		outfile = ".".join(tmpHat)
	else:
		# file is not gzipped
		tmp = (args.fin).split(".")
		n = len(tmp)
		tmpHat = tmp[0:(n-1)]
		tmpHat.append("umi")
		outfile = ".".join(tmpHat)

	rres = parse_umi(infile, outfile, ftype, gz)

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
# filters reads based on umi frequency. all distinct umis are counted and then 
# reads are written to 'outfile' that have UMIs which occur at least the 
# specified minimum number of times.
def parse_umi(infile, outfile, fmt, gz):

	if fmt == "fasta":
		rres = parse_umi_fasta(infile, outfile, gz)
	else:
		rres = parse_umi_fastq(infile, outfile, gz)

	return rres


def parse_umi_fasta(infile, outfile, gz):

	szout = ""
	
	if gz:
		p1 = sp.Popen("gunzip -c {}".format(infile).split(), stdout=sp.PIPE)
		fin = p1.stdout
	else:
		fin = open(infile, "r")
	
	num_read = 0
	num_written = 0
	
	for szl in fin:
		if szl[0] == ">":
			num_read += 1
			tmp = szl.strip().split(":")
			if len(tmp) > 1:
				num_written += 1
				szout += "{}\n".format(tmp[-1])
	
	fin.close()
	
	with open(outfile, "w") as fout:
		fout.write(szout)
	
	sys.stderr.write("read {} reads, wrote {} umi\n".format(num_read, num_written))
	
	return 0


def parse_umi_fastq(infile, outfile, min_count):

	idx = 0
	szout = ""
	
	if gz:
		p1 = sp.Popen("gunzip -c {}".format(infile).split(), stdout=sp.PIPE)
		fin = p1.stdout
	else:
		fin = open(infile, "r")
	
	num_read = 0
	num_written = 0
	
	for szl in fin:
		idx += 1
		if idx == 1:
			num_read += 1
			tmp = szl.strip().split(":")
			if len(tmp) > 1:
				num_written += 1
				szout += "{}\n".format(tmp[-1])
	
		elif idx > 3;
			idx = 0
	
	fin.close()
	
	with open(outfile, "w") as fout:
		fout.write(szout)
	
	sys.stderr.write("read {} reads, wrote {} umi\n".format(num_read, num_written))
	
	return 0



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

