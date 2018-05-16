#!/usr/bin/python
#==============================================================================
# bam-index-template.py
#
# Shawn Driscoll
# 20180118
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script is a framework for doing things with genome aligned reads. This
# code converts the alignments to SAM, if necessary, and indexes the reads.
#==============================================================================

import sys
import argparse
import math
import re
import os
import traceback
from collections import defaultdict
from time import localtime, time
from os.path import isfile
import hashlib
from random import random

#==============================================================================
# globals
#==============================================================================


#==============================================================================
# main
#==============================================================================

def main(args):

	is_pe = False
	is_bam = False
	tmp_base = ""
	aln_file = ""
	offset = 0
	read_index = None

	if not isfile(args.alignments):
		error_message("Input alignment file does not exists")
		return 1

	tmp_base = hashlib.md5(args.alignments).hexdigest()

	if re.search("\.bam$", args.alignments):
		is_bam = True
		message("Converting BAM alignments to SAM")
		aln_file = "{}.sam".format(tmp_base)
		rres = bam2sam(args.alignments, aln_file)

		if not isfile(aln_file):
			error_message("Conversion from BAM to SAM must have failed.")
			return 1

	else:
		aln_file = args.alignments


	##
	## open alignments and index reads
	##
	message("Indexing reads")
	read_index = build_read_index(aln_file)

	##
	## now we have the read names indexed and we can do stuff
	##
	for rname in read_index.keys():
		# do stuff
		pass


	if is_bam and isfile("{}.sam".format(tmp_base)):
		message("Removing temporary files")
		os.unlink("{}.sam".format(tmp_base))

	return 0

def build_read_index(sam_file):

	offset = 0
	read_index = defaultdict(list)

	with open(sam_file, "r") as fin:

		for szl in fin:
			if szl[0] != "@":
				# we have an alignment!
				aln = szl.strip().split("\t")
				qname = aln[0]
				if (int(aln[1]) & 0x1) != 0:
					is_pe = True

				read_index[qname].append(offset)

			offset += len(szl)

	message("{} distinct reads".format(len(read_index.keys())))

	return read_index


##
# function to call samtools for conversion of alignments
def bam2sam(infile, outfile):
	cmd = "samtools view -o {} {}".format(outfile, infile)
	rres = runcmd(cmd)
	return rres

def runcmd(cmd):
	sys.stderr.write("CMD: {}\n".format(cmd))
	return os.system(cmd)

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

parser = argparse.ArgumentParser(description="Template for processing alignments.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('alignments', type=str, 
	help="Alignments in SAM or BAM format for processeing.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

