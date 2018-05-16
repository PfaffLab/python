#!/usr/bin/env python
#==============================================================================
# gtf-get-contained.py
#
# Shawn Driscoll
# 20170111
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse a GTF to find transcripts and introns. Search for 
# features that are completely contained within a single intron of another feature.
# export the intron so that it can be added to a 'transcriptome' reference. 
# this will allow programs like kallisto to figure out if reads belong
# to the intron or the conatined feature. hopefully that may help eliminate
# the confusion between retained introns and these contained features.
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib
from multiprocessing import cpu_count

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

HBIN = 16000
STRAND_POS = "+"
STRAND_NEG = "-"

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	dgtf = {}
	dintrons = {}
	dtrange = {}
	lktable = {}
	idx = 0
	tid_wrapped = set()
	
	# load the GTF
	sys.stderr.write("Loading GTF {}\n".format(args.gtf))
	dgtf = parse_gtf(args.gtf)
	sys.stderr.write("Found {} transcripts\n".format(len(dgtf.keys())))
	
	# sort the exons within each transcript
	sys.stderr.write("Sorting exons within transcripts...\n")
	sort_parsed_gtf(dgtf)
	
	#
	# loop through and find all of the introns. also make a new set of regions
	# that are the full extent of each transcript. we'll also build a lookup
	# hash of the transcript regions
	#
	
	for tid in dgtf.keys():
		
		num_exons = len(dgtf[tid])
		
		# make a region for this transcript
		dtrange[tid] = Region(ref=dgtf[tid][0].rname, 
				start=dgtf[tid][0].start, end=dgtf[tid][num_exons-1].end, 
				tag=tid)
		dtrange[tid].set_id()
		
		# hash
		lhh = hash_region(dtrange[tid], HBIN)
		# check and add
		for hkey in lhh:
			if hkey not in lktable:
				lktable[hkey] = []
			# insert transcript id at this hash key
			lktable[hkey].append(tid)
		
		
		if num_exons > 1:
			# get introns
			
			for i in range(1, num_exons):
				r = Region(ref=dgtf[tid][i].rname, start=dgtf[tid][i-1].end+1, 
						end=dgtf[tid][i].start-1)
				r.set_id()
				
				if r.id not in dintrons:
					# add this intron
					dintrons[r.id] = r
					dintrons[r.id].trunk = []
				
				dintrons[r.id].trunk.append(tid)
	
	liid = dintrons.keys()
	
	#
	# now we loop through all of the introns and see if any of them contain 
	# a transcript. if it does then we need to gather them and determine
	# exactly how much of the intron to keep since usually the entire thing
	# will be overkill.
	bfound = False
	idx = 0
	for iid in dintrons.keys():
		# hash this bro
		lhh = hash_region(dintrons[iid], HBIN)
		bfound = False
		ileft = 1e9
		iright = -1
		# see if these hit any transcripts
		for hkey in lhh:
#			if bfound:
#				break
			
			if hkey in lktable:
				# check each transcript
				for tid in lktable[hkey]:
					if dintrons[iid].contains(dtrange[tid]):
						ileft = min([dtrange[tid].start-args.b, ileft])
						iright = max([dtrange[tid].end+args.b, iright])
						bfound = True
	
		if bfound:
			idx += 1
			# check ileft and iright boundaries. make sure they are inside of the intron
			ileft = max([ileft, dintrons[iid].start])
			iright = min([iright, dintrons[iid].end])
			
			if ileft < 1:
				ileft = 1

			dintrons[iid].start = ileft
			dintrons[iid].end = iright
			dintrons[iid].set_id()

			# print this intron out as a GTF row
			gout = GtfRow(rname=dintrons[iid].ref, start=dintrons[iid].start, 
						end=dintrons[iid].end, db="introns", type="exon")
			# get strand
			gout.strand = dgtf[dintrons[iid].trunk[0]][0].strand
			gout.gid = iid
			gout.tid = "INTID_{:08d}".format(idx)
			gout.attrs["gene_name"] = dgtf[dintrons[iid].trunk[0]][0].attrs["gene_name"]
			gout.attrs["oId"] = ",".join(dintrons[iid].trunk)
			for tmp in dintrons[iid].trunk:
				tid_wrapped.add(tmp)
			print "\t".join(map(str, gout.tolist()))
	
	# 
	# check for additional single-exon features that may not be within other gene loci
	#
	for tid in dgtf.keys():
		
		num_exons = len(dgtf[tid])
		
		if num_exons > 1:
			continue

		if tid in tid_wrapped:
			continue

		# we have another transcript to wrap
		idx += 1
		tid_wrapped.add(tid)
		ileft = max([1, dtrange[tid].start - args.b])
		iright = dtrange[tid].end + args.b
		gout = GtfRow(rname=dgtf[tid][0].rname, start=ileft, end=iright, 
			db="wrapper", type="exon")
		gout.strand = dgtf[tid][0].strand
		gout.gid = "{}:{}-{}".format(gout.rname, gout.start, gout.end)
		gout.tid = "INTID_{:08d}".format(idx)
		gout.attrs["oId"] = tid
		gout.attrs["gene_name"] = "{}_wrapper".format(dgtf[tid][0].attrs["gene_name"])

		print "\t".join(map(str, gout.tolist()))


	return 0


#==============================================================================
# classes
#==============================================================================

#
# class to hold a genomic region.
class Region(object):
	
	def __init__(self, ref="", start=0, end=0, strand="+", tag="", trunk=None):
		
		self.ref = ref
		self.start = int(start)
		self.end = int(end)
		self.strand = strand
		self.tag = tag
		# space for extra junk
		self.trunk = trunk
		self.id = ""
	
	def __str__(self):
		return "{}|{}:{}-{}|{}|{}".format(self.id, self.ref, self.start, self.end, self.strand, self.tag)

	#
	# set id
	def set_id(self):
		self.id = "{}:{}-{}".format(self.ref, self.start, self.end)
		return 0
	
	#
	# check for overlap between this region and another one
	def overlap(self, r):
				
		if r.ref != self.ref:
			return False
		
		# check for region overlap
		if r.start <= self.end and self.start <= r.end:
			return True
		
		return False
	
	#
	# check to see if this is contained by another
	def contained_by(self, r):
		
		if r.ref != self.ref:
			return False
		
		# check if r contains this
		if r.start <= self.start and r.end >= self.end:
			return True
	
		return False
	
	#
	# check if this region contains the other
	def contains(self, r):
		
		if r.ref != self.ref:
			return False
		
		# check if this contains r
		if self.start <= r.start and self.end >= r.end:
			return True
		
		return False

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
# defs
#==============================================================================

#
# parse a GTF into a dict of GtfRow objects
def parse_gtf(fname):
	# variables
	gtfdb = {}
	grow = None

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:

		grow = GtfRow()
		grow.parse(szl.strip())

		if grow.type != "exon":
			continue

		if grow.tid not in gtfdb:
			gtfdb[grow.tid] = []

		gtfdb[grow.tid].append(grow)

	fin.close()

	return gtfdb

def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	#
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

#
# function to sort the exons of a parsed gtf within each transcript
def sort_parsed_gtf(dgtf):
	
	idx = []
	lstarts = []
	
	for tid in dgtf.keys():
		idx = []
		lstarts = []
		
		for grow in dgtf[tid]:
			lstarts.append(grow.start)
		
		if len(lstarts) > 1:
			# check sort
			idx = np.argsort(lstarts)
			ltmp = list(dgtf[tid])
			for i in range(len(idx)):
				dgtf[tid][i] = ltmp[idx[i]]
	
	return 0


#--
# hash_region
# r is a list with [ref, start, end] for the region
# and bin is a binning integer for hashing.
# the region may fall into multiple bins. 
# returns a list of hashes
def hash_region(r, rbin):

	rbin = int(rbin)
	hstart = int(r.start)/rbin
	hend = int(r.end)/rbin

	hh = []
	for i in range(hstart, hend+1):
		hh.append("{}:{}".format(r.ref, i))

	return(hh)

#--
# region_overlap
# compares regions a and b to see if they overlap. returns 
# a list: [0/1, overlap length]
def region_overlap(a, b):
	astart = float(a[1])
	aend = float(a[2])
	bstart = float(b[1])
	bend = float(b[2])
	ovl = 0
	olen = 0

	# check overlap
	if aend >= bstart and astart <= bend:
		ovl = 1
		olen1 = aend-bstart
		olen2 = bend-astart
		len1 = aend-astart
		len2 = bend-bstart
		# overlap length is the minimum of all of these lengths
		olen = min([olen1, olen2, len1, len2])

	return([ovl, olen])

#--
# build_range_lookup_hash
# lpoints is a list of, at minimum, two element lists that shall be hashed
def build_range_lookup_hash(lranges, bin):

	rid = []
	bin = int(bin)
	didx = {}

	for i in range(len(lranges)):
		hh = hash_region(lranges[i], bin)

		for j in range(len(hh)):
			rid = hh[j]

			if rid not in didx:
				didx[rid] = []

			# insert into the hash
			didx[rid].append(list(lranges[i]))

	return(didx)


#--
# find_point_hash_hits
# look up a point or region in a point hash
def find_region_hash_hits(a, h, bin):

	# hash a, a region
	hh = hash_region(a, bin)
	# check for hits
	hids = {}

	for i in range(len(hh)):
		if hh[i] in h:
			# maybe, get list of elements at this node
			lcand = h[hh[i]]
			for j in range(len(lcand)):
				# overlap?
				rres = region_overlap(a, lcand[j])
				if rres[0]==1:
					# keep this hit
					hid = range_to_id(lcand[j])
					hids[hid] = 0

	return(hids.keys())



#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="...")
parser.add_argument('gtf', type=str, help="GTF annotation")
parser.add_argument('-b', type=int, default=200, action="store", 
	help="Padding around maximum length transcript found within an intron [200]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

