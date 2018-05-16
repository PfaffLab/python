#!/usr/bin/python
#==============================================================================
# build-gtf.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This builds GTF files from the annotation tables that one can download
# from UCSC. If there is a separate transcript to gene name map then 
# it can be specified with -a <filename> and it should be a two column (minimum)
# with the first column being transcript ids and the second being gene names
#  
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
	tid_to_gname = {}

	if args.a is not None:
		fin = open(args.a, "r")
		for szl in fin:
			aln = szl.strip().split("\t")
			tid_to_gname[aln[0]] = aln[1]
		fin.close()

	#
	# open input file
	fin = open(args.tsv, "r")
	flog = open("buildGtf_{}_droppedTooShort.txt".format(args.t), "w")

	dtid = {}
	for szl in fin:
		rres = tsv_row_to_gtf_rows(szl.strip(), args.t, tid_to_gname, dtid, args)
		if len(rres) > 0:
			for k in rres:
				print "\t".join(map(str, k))
		else:
			flog.write(szl)

	fin.close()
	flog.close()

	return 0

#
# function to reformat a single row from the input TSV into a full 
# transcript set of rows (minus last field) 
def tsv_row_to_gtf_rows(sz, dtype, ttg, dt, args):
	aln = sz.split("\t")
	lout = []
	tid_length = 0

	if args.kg:

		estarts = aln[8].split(",")
		eends = aln[9].split(",")

	else:

		estarts = aln[9].split(",")
		eends = aln[10].split(",")

	if args.kg:
		tid = aln[0]
	else:
		tid = aln[1]

	if tid not in dt:
		dt[tid] = 0
	else:
		dt[tid] += 1
		tid = "{}-{}".format(tid, dt[tid])
		
	# determine transcript length

	for i in range(len(estarts)):
		if len(estarts[i]) > 0 and len(eends[i]) > 0:
			tid_length += int(eends[i])-int(estarts[i])

	if tid_length > args.min_length:

		for i in range(len(estarts)):
			if len(estarts[i]) > 0:

				if args.kg:
					grow = GtfRow(rname=aln[1], db=dtype, type="exon", start=int(estarts[i])+1, end=int(eends[i]), strand=aln[2])
					grow.tid = tid
					grow.gid = aln[0]

					if aln[0] in ttg:
						grow.attrs["gene_name"] = re.sub("[\(\)\{\}\[\]\,\.\s]", "_", ttg[aln[0]])
					else:
						if len(aln[10]) > 0:
							grow.attrs["gene_name"] = re.sub("[\(\)\{\}\[\]\,\.\s]", "_", aln[10])
						else:
							grow.attrs['gene_name'] = "unknown"

				else:
					grow = GtfRow(rname=aln[2], db=dtype, type="exon", start=int(estarts[i])+1, end=int(eends[i]), strand=aln[3])
					grow.tid = tid
					grow.gid = aln[12]

					if aln[1] in ttg:
						grow.attrs["gene_name"] = re.sub("[\(\)\{\}\[\]\,\.\s]", "_", ttg[aln[1]])
					else:
						grow.attrs["gene_name"] = re.sub("[\(\)\{\}\[\]\,\.\s]", "_", aln[12])

				lout.append(grow.tolist())

	return lout


#==============================================================================
# GtfRow Class
#==============================================================================

class GtfRow(object):
	def __init__(self, rname="", db="", type="", start=0, end=0, strand="."):
		self.rname = rname
		self.db = db
		self.type = type
		self.start = start
		self.end = end
		self.strand = strand
		
		self.tid = ""
		self.gid = ""
		
		self.attrs = {}
		
	# parse info in GTF row, sz, into this object
	def parse(self, sz):
		aln = sz.strip().split("\t")
		self.rname = aln[0]
		self.db = aln[1]
		self.type = aln[2]
		self.start = int(aln[3])
		self.end = int(aln[4])
		self.strand = aln[6]
		
		# parse attributes from field 8
		fsplit = aln[8].split("\"")
		n = len(fsplit)-1
		i = 0
		while i < n:
			key = re.sub(';','',fsplit[i])
			self.attrs[key.strip()] = fsplit[i+1].strip()
			i += 2

		self.tid = self.attrs['transcript_id']
		self.gid = self.attrs['gene_id']
		
		return 0

	#--
	# tolist
	# return a list version of this with elements in place of the columns 
	# of a GTF
	def tolist(self):
		ltmp = [self.rname, self.db, self.type, 
			int(self.start), int(self.end), ".", self.strand, "."]
		
		szattr = "transcript_id \"{}\"; gene_id \"{}\";".format(self.tid, self.gid)
		akey = self.attrs.keys()
		if len(akey) > 0:
			# append additional attributes to the szattr string
			for aid in akey:
				szattr += " {} \"{}\";".format(aid, self.attrs[aid])
		
		ltmp.append(szattr)
		
		return(ltmp)


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('tsv', type=str, help="TSV file downloaded from UCSC used to build a GTF")
parser.add_argument('-a', type=str, default=None, 
	help="Transcript id to gene name table (usually Ensembl annotation has this)")
parser.add_argument('-t', type=str, default="gene", 
	help="Type field for GTF rows [gene]")
parser.add_argument('--kg', action="store_const", const=True, default=False, 
	help="Set this if the source table is 'knownGene' and the starts/ends of exons is in columns 9 and 10 instead of 10 and 11")
parser.add_argument('-l', '--min-length', type=int, default=50, action="store", 
	help="Minimum total transcript length [50]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

