#!/usr/bin/env python
#==============================================================================
# merge-gtf.py
#
# Shawn Driscoll
# 20160127
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script accepts multiple GTF files and combines them collapsing transcripts
# with identical intron chains into single transcripts. In those cases the 
# transcript with the longest overall genomic length is kept (i.e. they might
# have the same set of introns but different 3' or 5' UTR lengths).
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

#HOME = expanduser("~")

_CHROM = 0
_STRAND = 6
_START = 3
_END = 4
_TYPE = 2
_ATTR = 8

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	dtid = {}
	dintChain = {} # dict keyed by hashes of intron chains
	dsingleExons = {} # dict for single exon features

	dfinalTranscripts = {}

	# check the files
	for f in args.gtf_files:
		if not isfile(f):
			sys.stderr.write("Input file does not exist: {}".format(f))
			return 1

	# load the GTFs and establish which ones will be collapsed
	for f in args.gtf_files:
		sys.stderr.write("Loading {}\n".format(f))
		parse_gtf(f, dtid, dintChain, dsingleExons)


	sys.stderr.write("Building non-redundant transcript set...\n")

	# loop through multi-exon features

	for ichain in dintChain.keys():

		maxId = dintChain[ichain][0]
		maxLen = dtid[maxId]['length']

		if len(dintChain[ichain]) > 1:
			# multiple features, get one with longest sequence
			for tid in dintChain[ichain]:
				if dtid[tid]['length'] > maxLen:
					maxLen = dtid[tid]['length']
					maxId = tid

		chrom = dtid[maxId]['chrom']

		if chrom not in dfinalTranscripts:
			dfinalTranscripts[chrom] = []

		dfinalTranscripts[chrom].append([maxId, dtid[maxId]['start'], ",".join(dintChain[ichain])])

		#print dfinalTranscripts[chrom][-1]

	# loop through single-exon features
	
	for ichain in dsingleExons.keys():

		# if there is more than one transcript id in this bin they are all
		# the same length so we'll just take the first one
		maxId = dsingleExons[ichain][0]

		chrom = dtid[maxId]['chrom']

		if chrom not in dfinalTranscripts:
			dfinalTranscripts[chrom] = []

		dfinalTranscripts[chrom].append([maxId, dtid[maxId]['start'], ",".join(dsingleExons[ichain])])

	# loop through the chromosomes and produce the final GTF output

	# open output file
	fout = open(args.gtf_out, "w")

	for chrom in sorted(dfinalTranscripts.keys()):
		
		lchrom = dfinalTranscripts[chrom]

		# sort by starts
		lstarts = []
		for i in range(len(lchrom)):
			lstarts.append(lchrom[i][1])

		sort_order = np.argsort(lstarts)

		for j in range(len(sort_order)):
			i = sort_order[j]

			tid = lchrom[i][0]
			# print transcript rows
			for n in range(len(dtid[tid]['rows'])):
				lout = list(dtid[tid]['rows'][n])
				lout[_ATTR] += " oId \"{}\";".format(lchrom[i][2])
				fout.write("\t".join(lout) + "\n")

	fout.close()

	return 0

def parse_gtf(f, dtid0, dintChain0, dsingleExons0):

	lbuffer = []
	current_tid = ""
	last_tid = ""
	arl = []

	# open up the file and parse it out

	fin = open(f, "r")

	for szl in fin:
		arl = szl.strip().split("\t")

		# skip rows that are not exons
		if arl[_TYPE] != "exon":
			continue

		current_tid = get_tid(szl)

		if current_tid != last_tid:
			if len(lbuffer) > 0:
				dcurrent = parse_transcript(lbuffer)
				dtid0[last_tid] = dcurrent
				
				if not dcurrent['se']:
					# multi-exon feature
					if dcurrent['ichain'] not in dintChain0:
						dintChain0[dcurrent['ichain']] = []

					dintChain0[dcurrent['ichain']].append(last_tid)
				else:
					# single-exon feature
					if dcurrent['ichain'] not in dsingleExons0:
						dsingleExons0[dcurrent['ichain']] = []

					dsingleExons0[dcurrent['ichain']].append(last_tid)

			# new buffer
			lbuffer = []

		lbuffer.append(list(arl))
		last_tid = current_tid

	fin.close()

	# deal with last one

	if len(lbuffer) > 0:
		dcurrent = parse_transcript(lbuffer)
		dtid0[last_tid] = dcurrent

		if not dcurrent['se']:
			# multi-exon feature
			if dcurrent['ichain'] not in dintChain0:
				dintChain0[dcurrent['ichain']] = []

			dintChain0[dcurrent['ichain']].append(current_tid)
		else:
			# single-exon feature
			if dcurrent['ichain'] not in dsingleExons0:
				dsingleExons0[dcurrent['ichain']] = []

			dsingleExons0[dcurrent['ichain']].append(current_tid)

	return 0
		

def get_tid(sz):

	r = re.search("transcript_id \"([^\"]+)\"", sz)
	if r:
		return(r.group(1))

	return ""


#
# parse_transcript
# input: lrows - list of GTF rows associated with a single transcript
# this function will create the intron chain hash and sum the length of the 
# transcript
def parse_transcript(lrows):

	# variables
	dout = {}
	lstarts = []
	lmodrows = []
	lexons = []
	lintrons = []
	strand_pos = True
	ichain = ""
	tlen = 0

	n = len(lrows)

	# confirm strand
	if lrows[0][_STRAND] == "-":
		strand_pos = False

	# start the output dict
	dout['se'] = False
	dout['strand_pos'] = strand_pos

	# collect starts
	for i in range(n):
		lstarts.append(lrows[i][_START])

	# figure out sort order (order is from lowest to highest)
	sort_order = np.argsort(lstarts)

	# get transcript length
	for i in range(n):
		tlen += (int(lrows[i][_END])-int(lrows[i][_START])+1)

	ichain = "{}|".format(lrows[0][_CHROM])

	if n > 1:
		# build intron list
		for i in range(1, len(sort_order)):
			lintrons.append([int(lrows[sort_order[i-1]][_END])+1, int(lrows[sort_order[i]][_START])-1])

		for i in range(len(lintrons)):
			ichain += (":".join(map(str, lintrons[i])) + "|")

	else: 
		# single exon condition
		ichain += "{}:{}".format(lrows[0][_START], lrows[0][_END])
		dout['se'] = True

	# create final dict
	dout['ichain'] = ichain
	dout['length'] = tlen
	dout['rows'] = lrows
	dout['chrom'] = lrows[0][_CHROM]
	dout['start'] = int(lrows[sort_order[0]][_START])

	return(dout)




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

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('gtf_out', type=str, 
	help="Name of output GTF")
parser.add_argument('gtf_files', type=str, nargs='+', 
	help="Input file")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
#	except Exception as e:
#		sys.stderr.write("\nerror({}): {}\n".format(e.errno, e.strerror))

