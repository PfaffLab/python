#!/usr/bin/env python
#==============================================================================
# juncdb-parse-hits.py
#
# Shawn Driscoll
# 20160422
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Making a cleaner version of parse-junc-db-hits.py
#==============================================================================

import sys, argparse 
import math 
import re
import copy
from os.path import isfile, expanduser
import subprocess as sp
import os

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

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	szl = ""
	aln = []
	endid = []
	dhits = {}
	dgtf = {}

	# load the GTF
	sys.stderr.write("Loading the GTF {}\n".format(args.gtf))
	dgtf = parse_gtf(args.gtf)

	# parse the hits
	sys.stderr.write("Parsing hits from {}\n".format(args.fa))
	dhits = parse_alignments(args.fa, dgtf)

	sys.stderr.write("Building final table\n")

	print "\t".join(["chrom", "start", "end", "count", "id", "strand", "donor_ovl", "donor_use", "acc_ovl", "acc_use", 
		"psi5", "psi3", "theta5", "theta3", "gene_name", "gene_id", "transcript_id"])

	for tid in dgtf.keys():
		jpos = jid_to_pos(tid)
		endid = jid_to_daid(tid, dgtf[tid][0].strand)

		# set info
		jcount = JCount()
		jcount.chrom = jpos[0]
		jcount.start = int(jpos[1])
		jcount.end = int(jpos[2])
		jcount.strand = dgtf[tid][0].strand
		jcount.id = tid
		jcount.gname = dgtf[tid][0].attrs['gene_name']
		jcount.gid = dgtf[tid][0].gid
		jcount.tid = dgtf[tid][0].tid

		# set counts
		if tid in dhits:
			jcount.jcount = dhits[tid]['use']

		if endid[0] in dhits:
			jcount.donor_ovl = dhits[endid[0]]['ovl']
			jcount.donor_use = dhits[endid[0]]['use']

		if endid[1] in dhits:
			jcount.acc_ovl = dhits[endid[1]]['ovl']
			jcount.acc_use = dhits[endid[1]]['use']

		jcount.update()

		print "\t".join(map(str, jcount.tolist()))

	return 0


#==============================================================================
# classes
#==============================================================================

class JCount(object):
	def __init__(self):
		self.chrom = ""
		self.start = 0
		self.end = 0
		self.id = ""
		self.strand = "+"
		self.gname = ""
		self.gid = ""
		self.tid = ""
		self.jcount = 0
		self.donor_ovl = 0
		self.donor_use = 0
		self.acc_ovl = 0
		self.acc_use = 0
		self.psi5 = 0
		self.psi3 = 0
		self.theta5 = 0
		self.theta3 = 0

	def update(self):
		if self.donor_use > 0:
			self.psi5 = self.jcount*1.0/self.donor_use

		if self.donor_ovl > 0:
			self.theta5 = self.donor_use*1.0/(self.donor_use+self.donor_ovl)
		elif self.donor_use > 0:
			self.theta5 = 1

		if self.acc_use > 0:
			self.psi3 = self.jcount*1.0/self.acc_use
		
		if self.acc_ovl > 0:
			self.theta3 = self.acc_use*1.0/(self.acc_use+self.acc_ovl)
		elif self.acc_use > 0:
			self.theta3 = 1

		return 0

	def tolist(self):
		lout = [self.chrom, self.start, self.end, self.jcount, self.id, 
			self.strand, self.donor_ovl, self.donor_use, self.acc_ovl, self.acc_use, 
			self.psi5, self.psi3, self.theta5, self.theta3, self.gname, self.gid, self.tid]

		return lout

#
# class to hold GTF row objects. 
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
# functions
#==============================================================================

def jid_to_pos(jid):
	ltmp = jid.split(":")
	ltmp2 = ltmp[1].split("-")

	return([ltmp[0], ltmp2[0], ltmp2[1]])

def jid_to_daid(jid, strand):
	# explode the jid
	ltmp = jid.split(":")
	ltmp2 = ltmp[1].split("-")

	lid = "{}:{}".format(ltmp[0], ltmp2[0])
	rid = "{}:{}".format(ltmp[0], ltmp2[1])

	if strand == "-":
		lid += ":3p"
		rid += ":5p"
	else:
		# default to positive orientation
		lid += ":5p"
		rid += ":3p"

	return([lid, rid])

#
# parse a GTF into a dict of GtfRow objects
def parse_gtf(fname):
	# variables
	gtfdb = {}
	grow = None
	idx = 0

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:
		aln = szl.split("\t")

		if aln[2] != "exon" or aln[1] != "psi":
			continue

		grow = GtfRow()
		grow.parse(szl.strip())

		if grow.tid not in gtfdb:
			idx += 1
			gtfdb[grow.tid] = []

		gtfdb[grow.tid].append(grow)

	fin.close()

	sys.stderr.write("parsed {} transcripts\n".format(idx))
	
	return gtfdb

#
# parse_alignments
# this function parses the Seal hits and generates a dict of features hit with 
# counts
def parse_alignments(fname, dgtf):

	szl = ""
	ltargets = []
	lhits = []
	feature_hits = {}
	wfactor = 1
	num_tied = 0
	idx = 0

	# load the GTF
#	sys.stderr.write("Loading the GTF {}\n".format(args.gtf))
#	dgtf = parse_gtf(args.gtf)

	# parse the hits
	fin = open(args.fa, "r")
	for szl in fin:
		if not re.match("^>", szl):
			continue
		
		idx += 1
		if (idx % 1000000) == 0:
			sys.stderr.write("parsed {} aligned reads\n".format(idx))

		wfactor = 1.0
		ltargets, lhits = parse_hit_name(szl)
		
		num_tied = 1
		if len(ltargets) > 1:

			# take best only. hit counts are sorted in descending order already
			hmax = lhits[0]

			# increment through the hit counts to find how many share the sme
			# max match count
			i = 1
			while i < len(lhits):
				if lhits[i] < hmax:
					break

				i += 1

			# hit adjustment factor is 1/<number of best hit features>
			wfactor = 1.0/i
			num_tied = i

		# count the hits into the features hit adding their names to the dict as we go
		for j in range(0,num_tied):

			if ltargets[j] not in feature_hits:
				feature_hits[ltargets[j]] = {'use': 0, 'ovl': 0}

			# assign hit to the junction donor/acceptor use or ovl counts
			if is_theta(ltargets[j]):
				# target is a donor/acceptor overlap count
				feature_hits[ltargets[j]]['ovl'] += wfactor
			else:
				feature_hits[ltargets[j]]['use'] += wfactor
				
				# get strand oriented donor and acceptor ids
				endid = jid_to_daid(ltargets[j], dgtf[ltargets[j]][0].strand)
				
				# add them to the dict if not already added
				if endid[0] not in feature_hits:
					feature_hits[endid[0]] = {'use': 0, 'ovl': 0}
				
				if endid[1] not in feature_hits:
					feature_hits[endid[1]] = {'use': 0, 'ovl': 0}

				# increment their 'use' counts
				feature_hits[endid[0]]['use'] += wfactor
				feature_hits[endid[1]]['use'] += wfactor

	# final message
	sys.stderr.write("parsed {} aligned reads (done)\n".format(idx))
	
	fin.close()	

	return feature_hits

def is_theta(id):
	if re.search("[35]p$", id):
		return True

	return False

# ---
# parse a hit name from one of the FASTA lines output from Seal
def parse_hit_name(sz):

	aln = sz.strip().split("\t")

	dhits = {}
	ltargets = []
	lhits = []

	# just one feature means there's only a read name and no targets
	if len(aln)==1:
		return ltargets, lhits

	# loop through targets. because of how the index is made we'll have hits to many mers
	# of the same target so we need to bin them
	for i in range(1,len(aln)):
		alnx = aln[i].split("=")
		alnxName = alnx[1].split("|")

		for i in range(len(alnxName)):
			if alnxName[i] not in dhits:
				dhits[alnxName[i]] = 0

			dhits[alnxName[i]] += int(alnx[2])
	
	# pass to lists
	for tid in dhits.keys():
		ltargets.append(tid)
		lhits.append(dhits[tid])

	# sort by hit count
	if len(lhits) > 1:
		o = list(np.argsort(lhits))
		o.reverse()
		ltmpa = list(lhits)
		ltmpb = list(ltargets)
		for i in range(len(lhits)):
			lhits[i] = ltmpa[o[i]]
			ltargets[i] = ltmpb[o[i]]

	return ltargets, lhits



# --
# runcmd
# run a system level command in subprocess. optionally you can return the process.
# if the process isn't returned then the function waits for the process to finish
def runcmd(cmd, returnProcess=False):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Parse the output of seal for counts of kmers to targets.")
parser.add_argument('gtf', type=str, help="Juncdb GTF annotation")
parser.add_argument('fa', type=str, help="seal hits as FASTA")

args = parser.parse_args()

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

