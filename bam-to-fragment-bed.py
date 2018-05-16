#!/usr/bin/env python
#==============================================================================
# template.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script parses BAM alignments and print out a bed file with full-fragment
# size regions. Spliced alignments are skipped. 
#
# process:
# if pe only look at first mate (0x40). get fragment size (aln[8])
# if orientation is forward, get start by subtracing soft clipped bases
# from the position (aln[3] and cigar @ aln[5]). if orientation is reversed
# find start by iterating to the end of the alignment via the CIGAR 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

import subprocess as sp
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
	drefs = {}

	# open bam
	if re.search("\.bam$", args.bam):
		p1 = sp.Popen("samtools view -hF 0x4 {}".format(args.bam).split(), stdout=sp.PIPE)
	elif re.search("\.sam$", args.bam):
		p1 = sp.Popen("samtools view -hSF 0x4 {}".format(args.bam).split(), stdout=sp.PIPE)

	for szl in p1.stdout:

		if szl[0] == "@":
			# header line
			tmp = parse_sam_header_row(szl)

			if tmp["type"] == "SQ":
				# keep ref length so we can make sure reads don't extend out past the ends
				drefs[tmp["data"][0]] = int(tmp["data"][1])
			continue

		aln = parse_bam_alignment(szl)

		if aln_is_aln(aln):
			alr = aln_get_left_right(aln)
			cflags = aln_check_flags(aln)

			if cflags["second_mate"]:
				continue

			lout = [aln["rname"], alr[0], alr[1], "@" + aln["qname"]]
			
			# if the alignment is forward oriented then we only need the
			# left position and we'll add the fragment length

			flen = args.l
			if cflags["pe"]:
				# use the fragment length
				flen = int(aln["tlen"])
				if flen < 0:
					flen = flen*-1

			if cflags["reversed"]:
				# right side is fixed, extend to the left
				lout[1] = lout[2]-flen+1
				if lout[1] < 1:
					lout[1] = 1
			else:
				# extend to the right
				lout[2] = lout[1]+flen-1
				if lout[2] > drefs[aln["rname"]]:
					lout[2] = drefs[aln["rname"]]

			print "\t".join(map(str, lout))

	return 0



#==============================================================================
# defs
#==============================================================================


#--
# parse_sam_header_row
# parses a sam header line into a list
def parse_sam_header_row(sz):
	aln = sz.strip().split("\t")
	dh = {}
	x = aln[0]
	dh["type"] = x[1:len(x)]
	ll = []
	for i in range(1,len(aln)):
		tmp = aln[i].split(":")
		ll.append(tmp[1])

	dh["data"] = ll
	return(dh)


#--
# parse_bam_alignment
# parse a bam alignment into a dict
def parse_bam_alignment(sz):
	# vars
	aln = []
	fields = ["qname", "flag", "rname", "pos", "mapq", "cigar", 
		"rnext", "pnext", "tlen", "seq", "qual", "ext"]
	ext_idx = 11
	daln = {}

	# explode
	aln = sz.strip().split("\t")
	nfields = len(aln)

	for i in range(0, ext_idx):
		daln[fields[i]] = aln[i]

	daln["ext"] = list(aln[ext_idx:nfields])

	return(daln)


def aln_is_spliced(daln):
	rres = False
	rr = re.search("[0-9]+N", daln["cigar"])
	if rr:
		rres = True
	return rres

#--
# explode_cigar
# explode cigar string into operations and lengths
# note: uses re
def explode_cigar(cigar):
	
	# using findall we can blow this thing up all in one shot
	res = []
	if cigar != "*":
		res = re.findall("([0-9]+)([MIDNSHPX=]+)", cigar)
	
	return(res)

#--
# aln_get_left_right
# get the left and right most coordinates (adding soft clips back in)
# for the alignment. expect a dict input.
def aln_get_left_right(daln):

	res = [-1, -1]
	left0 = 0
	left = 0
	right = 0

	if daln['cigar'] != "*" and aln_is_aln(daln):
		left0 = int(daln['pos'])
		left = left0
		right = left0

		# adjust left for soft-clipping?
		rres = re.search("^([0-9]+)S", daln["cigar"])
		if rres:
			left -= int(rres.group(1))

		# get ops that change position
		rres = re.findall("([0-9]+)([MDN])", daln["cigar"])
		# find right
		for i in range(len(rres)):
			right += int(rres[i][0])

		# adjust right for soft-clipping?
		rres = re.search("([0-9]+)S$", daln["cigar"])
		if rres:
			right += int(rres.group(1))

		right -= 1

		res = [left, right]

	return res

#--
# get a dict of True/False for all of the possible flags in an 
# alignment.
def aln_check_flags(daln):
	rres = {}
	flag = int(daln["flag"])

	rres["pe"] = (flag & 0x1)
	rres["proper"] = (flag & 0x2)
	rres["unmapped"] = (flag & 0x4)
	rres["mate_unmapped"] = (flag & 0x8)
	rres["reversed"] = (flag & 0x10)
	rres["mate_reversed"] = (flag & 0x20)
	rres["first_mate"] = (flag & 0x40)
	rres["second_mate"] = (flag & 0x80)
	rres["secondary"] = (flag & 0x100)
	rres["filtered"] = (flag & 0x200)
	rres["duplicate"] = (flag & 0x400)
	rres["supaln"] = (flag & 0x800)

	return(rres)


#--
# aln_is_aln
# return true if the alignment is aligned
def aln_is_aln(daln):
	res = True
	if int(daln["flag"]) & 0x4:
		res = False
	return(res)


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('bam', type=str, help="BAM alignments")
parser.add_argument('-l', type=int, default=180, 
	help="Estimate fragment length (for single end. paired-end is inferred from alignments) [180]")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

