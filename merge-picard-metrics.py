#!/usr/bin/python
#==============================================================================
# merge-picard-metrics.py
#
# Shawn Driscoll
# 20170621
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Merges the results from multiple CollectMultiple and CollectRna blah 
# stats. This would match up with what I'm getting from 'parse-run picard-stats'
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
import traceback

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

ALN_MET = "alignment_summary_metrics"
RNA_MET = "rna_metrics"
QUAL_MET = "quality_by_cycle_metrics"

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	dsamples = build_sample_list(args.input)

	fout_mets = "combined_picard_stats.tsv"
	fout_cov = "combined_picard_rna_coverage.tsv"

	lstats_header = []
	lstats = []
	lcov_bin = []
	lcov = []
	lcov_header = ["norm_position"]

	for k in sorted(dsamples.keys()):
		if isfile(dsamples[k]['aln']) and isfile(dsamples[k]['rna']):
			rna_header, rna_stats, rna_cov = parse_rna_metrics(dsamples[k]['rna'])
			aln_header, aln_stats = parse_aln_metrics(dsamples[k]['aln'])
			if len(lstats_header)==0:
				lstats_header = ["SAMPLE"] + aln_header + rna_header

			lstats.append([k] + aln_stats + rna_stats)

			if len(lcov_bin)==0:
				lcov_bin = rna_cov[0]

			lcov.append(rna_cov[1])
			lcov_header.append(k)
	
	with open(fout_mets, "w") as fout:
		fout.write("\t".join(lstats_header) + "\n")
		for l in lstats:
			fout.write("\t".join(l) + "\n")

	with open(fout_cov, "w") as fout:
		fout.write("\t".join(lcov_header) + "\n")
		for i in range(len(lcov_bin)):
			lout = [lcov_bin[i]]
			for j in range(len(lcov)):
				lout.append(lcov[j][i])

			fout.write("\t".join(lout) + "\n")

	return 0

#
# parse out the metrics from the ALN_MET type
def parse_aln_metrics(fname):
	lheader = []
	lstats = []

	with open(fname, "r") as fin:
		while True:
			szl = fin.readline()
			if re.search("^## METRICS", szl):
				break

		szl = fin.readline()
		lheader = szl.strip().split("\t")
		szl = fin.readline()
		lstats = szl.strip().split("\t")
		if len(lstats) < len(lheader):
			# pad stats with NA
			while len(lstats) < len(lheader):
				lstats.append("NA")

	return lheader, lstats

def parse_rna_metrics(fname):
	lheader = []
	lstats = []
	lcoverage = [[], []]

	with open(fname, "r") as fin:
		while True:
			szl = fin.readline()
			if re.search("^## METRICS", szl):
				break


		szl = fin.readline()
		lheader = szl.strip().split("\t")
		szl = fin.readline()
		lstats = szl.strip().split("\t")

		if len(lstats) < len(lheader):
			# pad stats with NA
			while len(lstats) < len(lheader):
				lstats.append("NA")

		for szl in fin:
			if re.search("^normalized_position", szl):
				break

		for szl in fin:
			aln = szl.strip().split("\t")
			if len(aln)==2:
				lcoverage[0].append(aln[0])
				lcoverage[1].append(aln[1])

	return lheader, lstats, lcoverage

# 
# we need to merge the 
def build_sample_list(flist):
	# samples dict
	dsamples = {}

	for f in flist:
		tmp = f.split(".")
		basename = tmp[0]
		ext = tmp[1]

		if ext != ALN_MET and ext != RNA_MET and ext != QUAL_MET:
			continue

		if basename not in dsamples:
			dsamples[basename] = {"rna":"", "aln": "", "qual": ""}

		if ext == ALN_MET:
			dsamples[basename]['aln'] = f
		elif ext == RNA_MET:
			dsamples[basename]['rna'] = f
		elif ext == QUAL_MET:
			dsamples[basename]['qual'] = f

	return dsamples

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
	return 0


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.")
parser.add_argument('input', type=str, nargs="+", 
	help="Pass all metrics files regardless of what they are.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except Exception, e:
		print_exception()
		sys.exit(1)

