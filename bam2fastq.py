#!/usr/bin/python
#==============================================================================
# bam2fastq.py
#
# Shawn Driscoll
# 20180226
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# I'm writing this because the current reformatting programs are becoming 
# to picky. reformat.sh fails from some stupid MD tag error and samtools
# won't do anything unless the input file has a proper header. F that. 
#==============================================================================

import sys
import argparse
#import math
import re
from os.path import isfile, expanduser
#from collections import defaultdict
from time import localtime
from Basics import messages as ms
import subprocess as sp
import gzip

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

# sam alignment fields
SAM_QNAME = 0
SAM_FLAG = 1
SAM_RNAME = 2
SAM_POS = 3
SAM_MAPQ = 4
SAM_CIGAR = 5
SAM_RNEXT = 6
SAM_PNEXT = 7
SAM_TLEN = 8
SAM_SEQ = 9
SAM_QUAL = 10
SAM_FIRST_ATTR = 11

#==============================================================================
# main
#==============================================================================

def main(args):

	##	
	## check input files
	##
	
	if not isfile(args.bam):
		ms.error_message("Input file does not exist")
		return 1
	
	if args.o is not None:
		
		if isfile(args.o):
			ms.message("Output file exists. Overwriting.")
		
		if args.z:
			if not re.search("\.gz$", args.o):
				args.o += ".gz"
		
		else:
			if re.search("\.gz$", args.o):
				# auto enable gzip if the output file name has .gz at the end
				args.z = True
	
	rres = core(args)

	return 0


def core(args):

	##
	## open input file
	if re.search("\.bam$", args.bam):
		p = sp.Popen("samtools view {}".format(args.bam).split(), stdout=sp.PIPE)
		fin = p.stdout
	
	else:	
		fin = open(args.bam, "r")
	
	ocount = 0
	szout = ""
	
	if args.o is None:
		fout = sys.stdout
	else:
		if args.z:
			fout = gzip.open(args.o, "w")
		else:
			fout = open(args.o, "w")
	
	# do it
	for szl in fin:
		if szl[0] == "@":
			# skip header lines
			continue
		
		aln = szl.strip().split("\t")
		flag = int(aln[1])
		
		if flag & 0x100:
			# don't do anything with secondary alignments
			continue
		
		if flag & 0x4:
			if args.no_unal:
				# skip unaligned if requested
				continue
		else:
			if args.unal:
				# only requested unaligned so skip this one
				continue
		
		# now we can export. is the read flipped?
		rname = "@" + aln[SAM_QNAME]
		if flag & 0x10:
			seq = revcomp(aln[SAM_SEQ])
			tmp = list(aln[SAM_QUAL])
			tmp.reverse()
			qual = "".join(tmp) 
		else:
			seq = aln[SAM_SEQ]
			qual = aln[SAM_QUAL]
		
		szout += "\n".join([rname, seq, "+", qual])
		szout += "\n"
		
		ocount += 1
		if (ocount % 1000000) == 0:
			fout.write(szout)
			szout = ""
	
	fout.write(szout)
	
	fin.close()
	fout.close()
	
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

parser = argparse.ArgumentParser(description="Convert SAM format alignments to FASTQ", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('bam', type=str, help="Input SAM or BAM format file with or without header.")
parser.add_argument('-o', type=str, default=None, 
	help="Output file name. STDOUT otherwise.")

parser.add_argument('-z', action="store_const", const=True, default=False, 
	help="gzip compress output")

unal_group = parser.add_mutually_exclusive_group(required=False)
unal_group.add_argument('--no-unal', action="store_const", const=True, default=False, 
	help="Do not write reads flagged as unaligned.")
unal_group.add_argument('--unal', action="store_const", const=True, default=False, 
	help="Only extract unaligned reads.")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

