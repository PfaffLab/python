#!/usr/bin/env python
#==============================================================================
# sylamer-enrichment-counts.py
#
# Shawn Driscoll
# 20151113
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Wrapper for doing a kind of manual enrichment count. Sylamer claims that 
# redundant source sequences should be collapsed but it seems like that 
# should only be true of features with overlapping genomic ranges since 
# we usually work from a genome. This script accepts the genome FASTA name, 
# a bed file of all of the UTR boundaries and a subset file to test
# for enrichment by transcript names. 
#==============================================================================

import sys, argparse, hashlib, os

import subprocess as sp

# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
import rpy2.robjects as robjects
r = robjects.r

HOME="/home/colfax"
SYLAMER = HOME + "/opt/sylamer/source/MOFIR_final/Sylamer/sylamer-static-linux-x86"

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	wordSubset = len(args.w) > 0
	wordAnnot = len(args.a) > 0

	sylamerCounts = {}
	uwords = {}
	sout = ""
	serr = ""
	tmpFasta = hashlib.md5(args.genome).hexdigest() + ".fa"
	tmpBed = hashlib.md5(args.bed).hexdigest() + ".bed"
	tmpWords = ""

	if wordSubset:
		tmpWords = hashlib.md5(args.w).hexdigest() + ".txt"
		fin = open(args.w, "r")
		for szl in fin:
			aln = szl.strip().split("\t")
			uwords[aln[0]] = 0
		fin.close()
		fout = open(tmpWords, "w")
		for szl in uwords.keys():
			fout.write(szl + "\n")
		fout.close()

	if wordAnnot:
		wordInfo = {}
		fin = open(args.a, "r")
		for szl in fin:
			aln = szl.strip().split("\t")
			wordInfo[aln[0]] = list(aln[1:])
			wordInfoLen = len(aln)-1

		fin.close()

	subset = {}
	
	# step 0 - make a fasta file of the total

	# step 1 - get total counts from sylamer

	# step 2 - subset the bed file and generate a fasta of distinct regions

	# step 3 - get subsetted counts from sylamer

	# step 4 - put both results into a table so we can run stats in R

	# --
	# create FASTA for total counts
	# --

	rres = fasta_from_bed(args.bed, args.genome, tmpBed, tmpFasta)


	# --
	# call sylamer to count mers and parse them into a hash
	# --

	cmd = SYLAMER + " -fasta {} -k {} --none".format(tmpFasta, args.mer)
	if wordSubset:
		cmd = cmd + " -words {}".format(tmpWords)

	ferr = open("/dev/null", "w")
	p1 = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=ferr)

	for sout in p1.stdout:
		aln = sout.strip().split("\t")
		sylamerCounts[aln[1]] = [aln[5], aln[6]]

	p1.stdout.close()
	ferr.close()

	# remove the fasta file
	os.unlink(tmpFasta)

	# --
	# now create the subset
	# --

	# parse in the subset file then open the BED and print out lines matching transcripts
	# in the subset
	fin = open(args.subset, "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		subset[aln[0]] = 0
	fin.close()

	fin = open(args.bed, "r")
	fout = open("subset_" + tmpBed, "w")
	for szl in fin:
		aln = szl.strip().split("\t")
		if aln[3] in subset:
			fout.write(szl)

	fin.close()
	fout.close()

	rres = fasta_from_bed("subset_" + tmpBed, args.genome, tmpBed, tmpFasta)
	os.unlink("subset_" + tmpBed)

	# --
	# count mers in the subset
	# --

	cmd = SYLAMER + " -fasta {} -k {} --none".format(tmpFasta, args.mer)
	if wordSubset:
		cmd = cmd + " -words {}".format(tmpWords)

	ferr = open("/dev/null", "w")
	p1 = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=ferr)

	for sout in p1.stdout:
		aln = sout.strip().split("\t")
		sylamerCounts[aln[1]] += [aln[5], aln[6]]

	p1.stdout.close()
	ferr.close()

	# remove the fasta file
	os.unlink(tmpFasta)

	if wordSubset:
		os.unlink(tmpWords)

	rphyper = robjects.r['phyper']

	# --
	# print this result out
	# --

	lout = []
	lpvals = []

	for mer in sylamerCounts.keys():
		#print "\t".join([mer]+sylamerCounts[mer])
		merVals = map(int, sylamerCounts[mer])
		testVals = [merVals[2], merVals[0], merVals[1]-merVals[0], merVals[3]]
		res = rphyper(testVals[0], testVals[1], testVals[2], testVals[3], lower_tail=False)
		tmp = [mer]
		if wordAnnot:
			if mer in wordInfo:
				tmp += wordInfo[mer]
			else:
				tmp += ["na" for i in range(wordInfoLen)]

		tmp += sylamerCounts[mer]
		lpvals.append(res.r_repr())
		tmp.append(lpvals[-1])
		lout.append(list(tmp))


	# now create adjusted pvalues
	rpadjust = robjects.r['p.adjust']
	padj = rpadjust(robjects.FloatVector(lpvals), method="BH")

#	print "\t".join(["word", "word.count", "total.count", "subset.word.count", "subset.total.count", "pval", "padj"])
	for i in range(len(lout)):
		lout[i].append(padj[i])
		print "\t".join(map(str, lout[i]))

	return 0


def fasta_from_bed(sourceBed, sourceFasta, targetBed, targetFasta):

	# output file for collapsed bed file
	fout = open(targetBed, "w")
	# sort command
	cmd = "sort -k1,1 -k2,2n {}".format(sourceBed)
	psort = sp.Popen(cmd.split(), stdout=sp.PIPE)
	# bedtools command
	cmd = "bedtools merge -s -i stdin"
	pbedtools = sp.Popen(cmd.split(), stdin=psort.stdout, stdout=fout, stderr=sp.PIPE)
	psort.stdout.close()
	sout,serr = pbedtools.communicate()
	psort.wait()
	fout.close()

	# --
	# create fasta and create counts
	# --

	fout = open(targetFasta, "w")
	cmd = "bedtools getfasta -s -fi {} -bed {} -fo {}".format(sourceFasta, targetBed, targetFasta)
	p1 = sp.Popen(cmd.split(), stdout=fout, stderr=sp.PIPE)
	sout,serr = p1.communicate()
	fout.close()

	os.unlink(targetBed)

	return 0


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

# setup arguments

parser = argparse.ArgumentParser(description="Exports counts from sylamer for total and subset.")
parser.add_argument('genome', type=str, help="Genome FASTA source.")
parser.add_argument('bed', type=str, help="BED file of all sequences to extract for analysis.")
parser.add_argument('subset', type=str, help="Tab-delim file with first column containing transcript ids to subset the population for enrichment.")
parser.add_argument('mer', type=int, help="Kmer length to quantify.")
parser.add_argument('-w', type=str, default="",
	help="Tab-delim file with first column containing mers to check for enrichment.")
parser.add_argument('-a', type=str, default="", 
	help="Word annotation file with kmer words in the first column and any number of additional columns following.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

