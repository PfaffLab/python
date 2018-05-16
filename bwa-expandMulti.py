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

import sys, argparse, pysam, re

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



#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	read_length = 0
	paired = False
	unaligned_output_bam = ""
	expanded_output_bam = ""
	ref_hash = {}
	
	#----------- check input files

	if not file_exists(args.alignments):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.alignments)
		return 1

	#----------- set some variables

	tmp = args.alignments.split(".")
	
	expanded_output_bam = ".".join(tmp[0:(len(tmp)-1)] + ["expanded.bam"])
	
	if args.u:
		unaligned_output_bam = ".".join(tmp[0:(len(tmp)-1)] + ["unaligned.bam"])

	

	#----------- get to work

	if args.S:
		sf = pysam.Samfile(args.alignments,"r")
	else:
		sf = pysam.Samfile(args.alignments, "rb")

	# peak into the file to set read length and paired/or not paired
	
	i = 0
	for aln in sf:
		i += 1
		paired = aln.flag & 0x1
		read_length = len(aln.seq)
		
		if i == 1:
			break
	
	sf.close()

	# open file back up and open unaligned output file if requested
		
	if args.S:
		sf = pysam.Samfile(args.alignments,"r")
	else:
		sf = pysam.Samfile(args.alignments, "rb")
	
	if args.u:
		# open unaligned output bam
		sout_un = pysam.Samfile(unaligned_output_bam, "wb", template=sf)
	
	# open expanded multi-mappers bam
	sout_exp = pysam.Samfile(expanded_output_bam, "wb", template=sf)	
	
	# build hash of reference names
	i = 0 
	for rref in sf.references:
		ref_hash[rref] = i
		i += 1 
	
	# read through the data
	if paired:
		# deal with PE mappings
		sys.stderr.write("PE implementation not yet complete")
		pass
	else:
		# deal with SE mappings
		
		for aln in sf:
			
			# find XA tag - if it contains content then we work on it. otherwise if the
			# read is unmapped and args.u is True we write it to sout_un or if it is
			# mapped we write it to sout_exp
			
			if aln.flag & 0x4:
				# unaligned
				if args.u:
					sout_un.write(aln)
			
			else:
				# print original alignment
				sout_exp.write(aln)

				# check XA tag
				
				idx = get_tag_index(aln.tags, "XA")
				if idx >= 0:
					# XA tag is present
					
					# split the XA data
					xa_list = list(aln.tags[idx])[1].split(";")
					
					# need the forward compliment read
					seq = aln.seq
					qual = aln.qual
					if aln.flag & 0x10:
						qual_rev = qual
						tmp = rev_comp(seq)
						seq = tmp
						tmp = list(qual)
						tmp.reverse()
						qual = "".join(tmp)
					else:
						tmp = list(qual)
						tmp.reverse()
						qual_rev = "".join(tmp)
						
					# loop through xa list
					for maln in xa_list:
						laln = maln.split(",")
						if len(laln) == 4:
							# contains info - we'll get an empty single element at the end of these since
							# the last alignment in the XA tag will have a semi-colon at the end of it
							
							laln_cigar = transform_cigar(laln[2])
							
							# -- need to figure out if the cigar matches the read length
							cigar_rlen = cigar_read_length(laln_cigar)
							if cigar_rlen != read_length:
								sys.stderr.write("Warning: multi-map cigar read length does not match actual read length\n")
							
							tmp = list(laln[1])
							rev_aln = tmp[0] == "-"
							pos = int("".join(tmp[1:len(tmp)]))-1
							
							# build a new alignment
							a = pysam.AlignedRead()
							a.qname = aln.qname
							
							# set sequence
							if rev_aln:
								a.seq = rev_comp(seq)
								a.qual = qual_rev
							else:
								a.seq = seq
								a.qual = qual
							
							# set flag
							a.flag = 0x100
							if rev_aln:
								a.flag |= 0x10
							
							a.rname = ref_hash[laln[0]]
							a.pos = pos
							a.mapq = 0

							# deal with cigar
							a.cigar = laln_cigar
							
							a.mrnm = -1
							a.mpos = 0
							a.isize = 0
							
							a.tags = [("NM", int(laln[3]))]
							
							sout_exp.write(a)
												

	sf.close()
	sout_exp.close()
	if args.u:
		sout_un.close()

	return 0


# --
# transform string CIGAR to op, length list
def transform_cigar(sz_cig):
	op_translation = {'M':0, 'I':1, "D":2, "N":3, "S":4, "H":5, "P":6, "=":7, "X":8 }
	
	# split the cigar into ops and lengths
	lengths = re.split("[A-Z]", sz_cig)
	ops = re.split("[0-9]+", sz_cig)

	lengths = lengths[0:(len(lengths)-1)]
	ops = ops[1:]
	
	lcig = []
	
	for i in range(len(ops)):
		lcig += [(op_translation[ops[i]], int(lengths[i]))]
	
	return lcig

# calculate read length from cigar
def cigar_read_length(cig):	
	
	rlen = 0
	for op in cig:
		if op[0] == 0 or op[0] == 1 or op[0] == 4:
			rlen += op[1]
	
	return rlen 
	

def rev_comp(read):
	#
	# reverse compliment a string of ACTGN's
	#
	ttable = {"A":"T", "C":"G", "T":"A", "G":"C", "N":"N"}
	n = len(read)
	rev_read = ""
	while n > 0:
		n -= 1
		rev_read += ttable[read[n]]

	return rev_read

# --
# get the index of a specific tag. returns -1 if the tag doesn't exist.
def get_tag_index(tags, tname):

	result = [i for i, v in enumerate(tags) if v[0] == tname]
	if len(result) == 1:
		return result[0]

	return -1

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
parser.add_argument('alignments', type=str, help="BAM or SAM alignments (use -S if SAM)")
parser.add_argument('-S', dest="S", action="store_const", const=True, default=False, help="Input is SAM (default: bam)")
parser.add_argument('--no-discordant', dest="no_discordant", action="store_const", const=True, default=False, 
				help="Do not expand or report paired alignments to different references (default: off)")
parser.add_argument('--no-singletons', dest="no_singletons", action="store_const", const=True, default=False,
				help="Do not expand singleton mappings (for PE data) (default: off)")
parser.add_argument('-u', dest="u", action="store_const", const=True, 
				default=False, help="Write unmapped reads to BAM output. (default: False)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
