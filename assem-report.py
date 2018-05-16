#!/usr/bin/python
#==============================================================================
# assem-report.py
#
# Shawn Driscoll
# 20170418
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Assembly report from Cufflinks and Stringtie attempts at assembling 
# alignments into transcripts compared to simulated expected
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

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

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	codeCounts = {}
	novel = []
	numq = 0
	expected = False
	found = False
	rres = []

	# load the control counts
	dexpr, dexpr0 = load_counts(args.counts, 5)

	# open cuffcmp tracking file
	fin = open(args.tracking, "r")
	for szl in fin:
		expected = False
		found = False
		aln = szl.strip().split("\t")

		if numq==0:
			numq = len(aln)-4
			for i in range(numq):
				# lists for discovered transcripts
				rres.append(Report())

		if aln[2] != "-":

			lhit = aln[2].split("|")
			tid = lhit[1]

			if tid in dexpr:
				expected = True
		else:
			tid = "novel"

		ccode = aln[3]

		# see what the assemblers had to say about this

		hits = []
		for i in range(4, len(aln)):
			r = parse_hit(aln[i])
			if r > 0:
				found = True
			hits.append(r)


		for i in range(numq):
			if hits[i] > 0:
				# found in this one

				if ccode != rres[i].codeCounts:
					rres[i].codeCounts[ccode] = 0
					rres[i].codeExpr[ccode] = 0

				rres[i].codeCounts[ccode] += 1
				rres[i].codeExpr[ccode] += hits[i]

				if expected:

					if ccode=="=":
						# good!
						rres[i].tp += 1
						rres[i].tpExpr += hits[i]
						rres[i].tidMatched.update([tid])

					else:
						# expected transcript but not a recovered assembly
						if ccode not in rres[i].wrongAssem:
							rres[i].wrongAssem[ccode] = 0
							rres[i].wrongAssemExpr[ccode] = 0

						rres[i].wrongAssem[ccode] += 1
						rres[i].wrongAssemExpr[ccode] += hits[i]

				else:
					# not expected
					if tid == "novel":
						rres[i].novel += 1
						rres[i].novelExpr += hits[i]

					else:
						# normal transcript id but it shouldn't have been detected
						if ccode=="=":
							# perfect match!
							rres[i].tidFalse.update([tid])

						if ccode not in rres[i].fp:
							rres[i].fp[ccode] = 0
							rres[i].fpExpr[ccode] = 0

						rres[i].fp[ccode] += 1
						rres[i].fpExpr[ccode] += hits[i]


	fin.close()

	# eval negatives
	totalExpected = len(dexpr.keys())
	totalNotExpected = len(dexpr0.keys())

	for i in range(numq):
		rres[i].totalPos = totalExpected
		rres[i].totalNeg = totalNotExpected

		print rres[i]


	return 0

# get either the expression or a -1 value if not assembled
def parse_hit(sz):
	rres = -1

	if sz != "-":
		aa = sz.split("|")
		rres = float(aa[3])

	return rres


def load_counts(f, min_count):
	dexpr = {} # for expressed transcripts
	dexpr0 = {} # for non-expressed

	fin = open(f, "r")
	# skip header
	szl = fin.readline()
	for szl in fin:
		aln = szl.strip().split("\t")
		if float(aln[3]) < min_count:
			dexpr0[aln[0]] = 0
		else:
			dexpr[aln[0]] = [float(aln[3]), float(aln[4]), 0]

	fin.close()
	return dexpr, dexpr0


#
# class to hold Report
class Report(object):
	def __init__(self):

		self.codeCounts = {}
		self.codeExpr = {}
		# expected transcript ids. some match annotation assembly some do not
		self.tp = 0
		self.tpExpr = 0
		self.wrongAssem = {}
		self.wrongAssemExpr = {}
		# not expected transcript ids..
		self.fp = {}
		self.fpExpr = {}
		self.tn = 0
		self.fn = 0
		self.novel = 0
		self.novelExpr = 0

		self.tidMatched = set()
		self.tidFalse = set()

		self.totalPos = 0
		self.totalNeg = 0

	
	def __str__(self):
		lout = self.tolist()
		return "\t".join(map(str, lout))

	def tolist(self):
		totalFp = 0
		totalFpExpr = 0
		fpString = ""
		fpExprString = ""

		totalWrongAssem = 0
		totalWrongAssemExpr = 0
		wrongAssemString = ""
		wrongAssemExprString = ""

		exactFp = 0
		if "=" in self.fp:
			exactFp = self.fp["="]

		for k in self.fp.keys():
			totalFp += self.fp[k]
			totalFpExpr += self.fpExpr[k]
			fpString += "{}:{:d}|".format(k, self.fp[k])
			fpExprString += "{}:{:0.4f}|".format(k, self.fpExpr[k])

		for k in self.wrongAssem.keys():
			totalWrongAssem += self.wrongAssem[k]
			totalWrongAssemExpr += self.wrongAssemExpr[k]
			wrongAssemString += "{}:{:d}|".format(k, self.wrongAssem[k])
			wrongAssemExprString += "{}:{:0.4f}|".format(k, self.wrongAssemExpr[k])

		prec1 = self.tp*1.0/(self.tp+exactFp)
		prec2 = self.tp*1.0/(self.tp+totalFp)
		recall = self.tp*1.0/(self.tp+self.fn)

#		lout = [fmt_float(prec1), fmt_float(prec2), fmt_float(recall), 
#			self.tp, exactFp, totalFp, self.tn, self.fn, self.tpExpr, totalFpExpr]


		# strict precision. positive is perfect assembly.
		tp = self.tp
		fp = self.fp["="]

		prec = tp*1.0/(tp+fp)
		recall = tp*1.0/self.totalPos

		# ratio of expected that were detected and correctly assembled
		prec2 = tp*1.0/(tp+totalWrongAssem)

		#lout = [self.tp, exactFp, totalWrongAssem, totalFp, self.tn, self.fn, self.novel]
		lout = [tp, fp, fmt_float(prec), fmt_float(recall), fmt_float(prec2)]

		# strict precision. positive is perfect assembly.
		tp = self.tp + totalWrongAssem
		fp = totalFp

		prec = tp*1.0/(tp+fp)
		recall = tp*1.0/self.totalPos

		lout += [tp, fp, fmt_float(prec), fmt_float(recall)]

		return lout


def fmt_float(v):
	return "{:0.4f}".format(v)

def collapse_keys_vals(d):

	lout = []
	for k in sorted(d.keys()):
		lout.append("{}:{}".format(k, d[k]))

	return "|".join(lout)




#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="..")
parser.add_argument('tracking', type=str, help="Cuffcompare 'cuffcmp.tracking' output.")
parser.add_argument('counts', type=str, help="Transcript read-counts from simulated counts.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

