#!/usr/bin/env python
#==============================================================================
# bwa-expandMulti.py
#
# Shawn Driscoll
# 20140113
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Expands multi-alignments from bwa 'aln' 'samse/sampe' output for use with
# tools like eXpress and RSEM.
#==============================================================================

import sys, argparse, re
import subprocess as sp

import pysam # remove after done with new main function


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

# -- sam fields
QNAME = 0
FLAG = 1
RNAME = 2
POS = 3
MAPQ = 4
CIGAR = 5
RNEXT = 6
PNEXT = 7
ISIZE = 8
SEQ = 9
QUAL = 10
FIRST_EXT = 11 

#==============================================================================
# main
#==============================================================================

def main(args):
	
	# -- variables
	
	rlen = 0
	paired = False
	unaligned_output_file = ""
	expanded_output_file = ""
	streaming = args.input == "-"
	iline = 0
	
	# -- set variables
	
	if streaming:
		unaligned_output_file = "bwa_unaligned"
		expanded_output_file = "bwa_expanded"
	else:
		tmp = args.input.split(".")
		unaligned_output_file = ".".join(tmp[0:(len(tmp)-1)] + ["unaligned"])
		expanded_output_file = ".".join(tmp[0:(len(tmp)-1)] + ["expanded"])
	
	if args.s:
		# output will be SAM
		unaligned_output_file += ".sam"
		expanded_output_file += ".sam"
	else:
		unaligned_output_file += ".bam"
		expanded_output_file += ".bam"
	
	# -- open files
	
	if streaming:
		fin = sys.stdin
	else:
		# if input is SAM just open it. if input is BAM we need to make a pipe
		if args.S:
			fin = open(args.input, "r")
		else:
			p1 = sp.Popen("samtools view -h {:s}".format(args.input).split(), stdout=sp.PIPE)
			fin = p1.stdout
	
	if args.s:
		# writing SAM output so we just need standard output streams
		sout_exp = file(expanded_output_file, "w")
		if args.u:
			sout_un = file(unaligned_output_file, "w")
	
	else:
		# writing BAM files so we have to setup a pipe
		p2 = sp.Popen("samtools view -bS -o {:s} -".format(expanded_output_file).split(), stdin=sp.PIPE)
		sout_exp = p2.stdin
		
		if args.u:
			p3 = sp.Popen("samtools view -bS -o {:s} -".format(unaligned_output_file).split(), stdin=sp.PIPE)
			sout_un = p3.stdin
	
	
	# read through input
	for szl in fin:
		
		if szl[0] == "@":
			# header line, write to output(s)
			sout_exp.write(szl)
			if args.u:
				sout_un.write(szl)
			
		else:
			iline += 1
			
			aln = szl.strip().split("\t")
			
			paired = int(aln[FLAG]) & 0x1
			rlen = len(aln[SEQ])
			
			if paired:
				pass
			
			else:
			
				# alignment row - is it aligned?
				if int(aln[FLAG]) & 0x4:
					# unaligned - write to unaligned file or not
					if args.u:
						sout_un.write(szl)
				
				else:
					# aligned read, write to main file
					# sout_exp.write(szl)
					# get the NM tag value
					# nm_val = get_tag_value(aln, "NM")
					# aln2 = list(aln[0:FIRST_EXT]) + ["NM:i:" + nm_val]
					# sout_exp.write("\t".join(aln2) + "\n")
					aln2 = strip_tag(list(aln), "XA")
					sout_exp.write("\t".join(aln2) + "\n")
					
					# check for XA tag
					xa_value = get_tag_value(aln, "XA")
					
					if len(xa_value) > 0:
						# has a value
						xa_list = xa_value.split(";")
						
						# need the forward compliment read
						seq = aln[SEQ]
						qual = aln[QUAL]
						
						if int(aln[FLAG]) & 0x10:
							# base alignment is reversed
							qual = rev_qual(qual)
							seq = rev_comp(seq)
							
						# loop through xa list
						for maln in xa_list:
							laln = maln.split(",")
							
							if len(laln) == 4:
								# contains info - we'll get an empty single element at the end of these since
								# the last alignment in the XA tag will have a semi-colon at the end of it
																
								# -- need to figure out if the cigar matches the read length
								cigar_rlen = cigar_read_length(laln[2])
								if cigar_rlen != rlen:
									sys.stderr.write("Warning: multi-map cigar read length does not match actual read length (line: {:d})\n".format(iline))
								
								# figure out the orientation of the multi-map location
								tmp = list(laln[1])
								rev_aln = tmp[0] == "-"
								pos = int("".join(tmp[1:len(tmp)]))
								
								# build a new alignment from original one
								
								aln_new = list(aln[0:(FIRST_EXT)]) + ["NM:i:" + laln[3]]
								aln_new[FLAG] = 0x100
								if rev_aln:
									aln_new[FLAG] |= 0x10
									
								aln_new[RNAME] = laln[0]
								
								if rev_aln:
									aln_new[SEQ] = rev_comp(seq)
									aln_new[QUAL] = rev_qual(qual)
								else:
									aln_new[SEQ] = seq
									aln_new[QUAL] = qual
								
								aln_new[POS] = str(pos)
								aln_new[MAPQ] = 0
								aln_new[CIGAR] = laln[2]
								aln_new[PNEXT] = 1
									
								sout_exp.write("\t".join(map(str, aln_new)) + "\n")

	
	# when done we have to close everything
	fin.close()
	sout_exp.close()
	if args.u:
		sout_un.close()
	
	return 0


def parse_cigar(sz_cig):
	# split the cigar into ops and lengths
	lengths = re.split("[A-Z]", sz_cig)
	ops = re.split("[0-9]+", sz_cig)

	lengths = lengths[0:(len(lengths)-1)]
	ops = ops[1:]

	return (ops, map(int, lengths))

# calculate read length from cigar
def cigar_read_length(sz_cig):	
	
	rlen = 0

	ops, lengths = parse_cigar(sz_cig)
	
	for i in range(len(ops)):
		if ops[i] == "M" or ops[i] == "S" or ops[i] == "I":
			rlen += lengths[i]
	
	return rlen 


def rev_qual(qual):
	tmp = list(qual)
	tmp.reverse()
	return "".join(tmp)

def rev_comp(read):
	#
	# reverse compliment a read
	#

	r2 = list(read.upper())
	r2.reverse()
	n = len(r2)
	
	for i in range(len(r2)):
		if r2[i] == "A":
			r2[i] = "T"
		elif r2[i] == "T":
			r2[i] = "A"
		elif r2[i] == "C":
			r2[i] = "G"
		elif r2[i] == "G":
			r2[i] = "C"
	
	return "".join(r2)

# --
# Return the value from an extended tag by tag name. If the tag isn't
# found then an empty string is returned.
def get_tag_value(aln, tname):

	n = len(aln)
	for i in range(FIRST_EXT, n):
		tag = aln[i].split(":")
		if tag[0] == tname:
			return(tag[2])
	
	return "" 

def strip_tag(aln, tname):

	n = len(aln)
	for i in range(FIRST_EXT, n):
		tag = aln[i].split(":")
		if tag[0] == tname:
			break
	
	aln.pop(i)
	
	return aln
	

# -- 
# checks that a file exists on disk. returns False if not.
def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True


#
# alignment accepted
#
def alignment_accepted(aln, args):
	
	result = True
	
	
	
	return result

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Expand multi-mappings from XA tag in BWA's samse/sampe output.")
parser.add_argument('input', type=str, help="Input reads (SAM, with -S, or BAM). Streaming SAM format is accepted as -")
parser.add_argument('-S', dest="S", action="store_const", const=True, default=False, help="Input is SAM (default: BAM)")
parser.add_argument('-s', dest='s', action="store_const", const=True, default=False, help="Ouptut SAM format (default: BAM)")
parser.add_argument('--no-discordant', dest="no_discordant", action="store_const", const=True, default=False, 
				help="Do not expand or report paired alignments to different references (default: off)")
parser.add_argument('--no-singletons', dest="no_singletons", action="store_const", const=True, default=False,
				help="Do not expand singleton mappings (for PE data) (default: off)")
parser.add_argument('-u', dest="u", action="store_const", const=True, 
				default=False, help="Write unmapped reads to BAM output. (default: False)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
