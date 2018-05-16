#!/usr/bin/env python
#==============================================================================
# gtf-to-chip-annotation.py
#
# Shawn Driscoll
# 20151022
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse gtf and create a bed format file containing positions of TSS, 
# boundaries of gene loci and boundaries of promoters
#==============================================================================

import sys, argparse, math, re
import numpy as np

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables

	# parse GTF into lookup hashes
	sys.stderr.write("Parsing GTF...\n")
	ltranscripts = load_gtf(args.gtf_file, args.p_bound_up, args.p_bound_down)

	for i in range(len(ltranscripts)):
		for j in range(len(ltranscripts[i])):
			print "\t".join(map(str, ltranscripts[i][j]))


	return 0

# --
# function to parse the GTF and create hashes of all features
def load_gtf(f, pboundUp, pboundDown):

	# list of parsed transcripts (each element is a list that includes TSS, promoter region, exons, introns and overall transcript boundary)
	ltranscripts = []
	lbuffer = []
	current_tid = ""
	last_tid = ""

	# open up the file and parse it out

	fin = open(f, "r")

	for szl in fin:
		arl = szl.strip().split("\t")

		current_tid = get_tid(szl)

		if current_tid != last_tid:
			if len(lbuffer) > 0:
				ltranscripts.append(parse_transcript(lbuffer, pboundUp, pboundDown))

			# new buffer
			lbuffer = []

		lbuffer.append(list(arl))
		last_tid = current_tid

	fin.close()

	# deal with last one

	if len(lbuffer) > 0:
		ltranscripts.append(parse_transcript(lbuffer, pboundUp, pboundDown))

	return(ltranscripts)
		

def get_tid(sz):

	r = re.search("transcript_id \"([^\"]+)\"", sz)
	if r:
		return(r.group(1))

	return ""

	
#
# parse_transcript
# input: lrows - list of GTF rows associated with a single transcript
# output: list of properties and boundaries in a bed-like format
def parse_transcript(lrows, pup, pdown):

	# variables
	lout = []
	lmodrows = []
	lstarts = []
	lexons = []
	lintrons = []
	ltss = []
	lpromoter = []
	strand_pos = True

	# boundaries for promotor region
	promotor_bound_pos = [-pup, pdown]
	promotor_bound_neg = [-pdown, pup]

	n = len(lrows)

	# reformat rows
	for i in range(n):
		lmodrows.append(reformat_gtf_row(lrows[i]))

	# confirm strand
	if lmodrows[0]['strand'] == "-":
		strand_pos = False

	# collect starts
	for i in range(n):
		lstarts.append(lmodrows[i]['start'])

	# figure out sort order (order is from lowest to highest)
	sort_order = np.argsort(lstarts)

	# build exon list
	for i in sort_order:
		lexons.append([lmodrows[i]['start'], lmodrows[i]['end']])

	if len(lexons) > 1:
		# build intron list
		for i in range(1, len(sort_order)):
			lintrons.append([lmodrows[sort_order[i-1]]['end']+1, lmodrows[sort_order[i]]['start']-1])

	if strand_pos:
		# positive strand so tss is the start of the first exon and promoter will be based on that
		ltss.append([lmodrows[sort_order[0]]['start']])
		lpromoter.append([
			lmodrows[sort_order[0]]['start']+promotor_bound_pos[0], 
			lmodrows[sort_order[0]]['start']+promotor_bound_pos[1] ])
	else:
		# negative strand so tss is the end of the last exon and promoter will be based on that
		ltss.append([lmodrows[sort_order[-1]]['end']])
		lpromoter.append([
			lmodrows[sort_order[-1]]['end']+promotor_bound_neg[0], 
			lmodrows[sort_order[-1]]['end']+promotor_bound_neg[1] ])


	# get ready to make a single list out of this thing

	# get tid, gid and gene name if present
	tid = lmodrows[0]['transcript_id']
	gid = lmodrows[0]['gene_id']
	gname = gid
	strand = lmodrows[0]['strand']
	chrom = lmodrows[0]['chrom']

	if "gene_name" in lmodrows[0]:
		gname = lmodrows[0]["gene_name"]

	# first the transcript boundary
	lout.append([chrom, 
		lmodrows[sort_order[0]]['start']-1, 
		lmodrows[sort_order[-1]]['end'], 
		tid, 0, strand, 
		gid, gname, "transcript"])

	# now the tss
	lout.append([chrom, 
		ltss[0][0]-1, 
		ltss[0][0], 
		tid, 0, strand, gid, gname, "tss"])

	# now the promoter
	lout.append([chrom, 
		lpromoter[0][0]-1, lpromoter[0][1], 
		tid, 0, strand, gid, gname, "promoter"])

	# now add the exon(s)
	
	if strand_pos:
		enum = 1
	else:
		enum = len(lexons)

	for i in range(len(lexons)):
		lout.append([chrom, 
			lexons[i][0]-1, lexons[i][1], 
			tid, 0, strand, gid, gname, "exon_{}_of_{}".format(enum, len(lexons)) ])

		if strand_pos:
			enum += 1
		else:
			enum -= 1


	# now add the introns(s)
	
	if len(lintrons) > 0:

		if strand_pos:
			enum = 1
		else:
			enum = len(lintrons)

		for i in range(len(lintrons)):
			lout.append([chrom, 
				lintrons[i][0]-1, lintrons[i][1], 
				tid, 0, strand, gid, gname, "intron_{}_of_{}".format(enum, len(lintrons)) ])

			if strand_pos:
				enum += 1
			else:
				enum -= 1

	# make sure there are no start positions less than 0

	for i in range(len(lout)):
		if lout[i][1] < 0:
			lout[i][1] = 0

	return lout

#
# reformats a tab exploded list gtf row into a dict so we can grab stuff more easily
def reformat_gtf_row(lrow):

	# start with the attributes field
	drow = parse_gtf_attr(lrow[8])

	drow['chrom'] = lrow[0]
	drow['start'] = int(lrow[3])
	drow['end'] = int(lrow[4])
	drow['strand'] = lrow[6]
	drow['feature'] = lrow[2]

	return(drow)


# --
# parse the attributes field of a GTF row into a hash
def parse_gtf_attr(sz):
	# split string on ;
	s1 = sz.split(";")
	attr = {}

	for i in range(len(s1)):
		s1[i] = s1[i].strip()
		if len(s1[i]) > 0:
			m = re.search("^([^\s]+) \"([^\"]+)\"", s1[i])
			attr[m.group(1)] = m.group(2)

	return attr
		

def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('gtf_file', type=str, help="Input GTF file")
parser.add_argument("--p-bound-up", type=int, default=1000, action="store", 
	help="Upstream distance from TSS for promoter boundary [1000]")
parser.add_argument("--p-bound-down", type=int, default=100, action="store", 
	help="Downstream distance from TSS for promoter boundary [100]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

