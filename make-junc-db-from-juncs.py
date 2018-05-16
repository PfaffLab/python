#!/usr/bin/env python
#==============================================================================
# make-junc-db-from-juncs.py
#
# Shawn Driscoll
# 20160415
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This version of make-junc-db accepts junctions such as those produced
# by bam2junc or juncs-and-boundaries. A GTF may be provided to annotate
# and include known junctions.
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
	
	juncs = []
	djuncs = {}
	dda = {} # donor/acceptor dict
	grow = None
	ngid = ""
	ngidx = 0
	juncdb_gtf = ""

	#	
	# set some values
	#

	juncdb_gtf = "{}.{}.gtf".format(args.stub, args.l)
	juncdb_fa = "{}.{}.fa".format(args.stub, args.l)
	juncdb_index = "{}.{}.index.fa".format(args.stub, args.l)

	# load annotation if provided. 
	if args.gtf is not None:
		sys.stderr.write("Loading junctions from {}\n".format(args.gtf))
		juncs = parse_gtf_to_junctions(args.gtf)
		
	# load junctions from file
	sys.stderr.write("Loading junctions from {}\n".format(args.juncs[0]))
	juncs += parse_juncs(args.juncs[0], args.s)
	if len(args.juncs) > 1:
		for i in range(1, len(args.juncs)):
			sys.stderr.write("Loading junctions from {}\n".format(args.juncs[i]))
			juncs += parse_juncs(args.juncs[i], args.s)
	
	
	# build a dict from the loaded junctions
	sys.stderr.write("Merging junctions from all inputs\n")
	djuncs = junc_list_to_dict(juncs)

	# if user provided a gtf we can use the 5' and 3' ends of annotated junctions
	# to attempt to annotate novel junctions
	
	if args.gtf is not None:
		sys.stderr.write("Annotating novel junctions by propagating donor/acceptor annotation\n")
		# build a dict with 5' and 3' ends of junctions from only those with annotation
		for jid in djuncs.keys():
			# is this one annotated?
			junc = djuncs[jid]
			if len(junc.tid) > 0:
				# yes - insert a key
				did = junc.get_donor_id()
				if did not in dda:
					dda[did] = [junc.tid, junc.gid, junc.gname, junc.strand]
					
				aid = junc.get_acc_id()
				if aid not in dda:
					dda[aid] = [junc.tid, junc.gid, junc.gname, junc.strand]
		
		# now flip through the dict again and check for unannotated junctions to 
		# see if we can annotate them
		for jid in djuncs.keys():
			junc = djuncs[jid]
			if len(junc.tid)==0:
				# let's check
				ld = []
				la = []
				
				# check for donor
				did = junc.get_donor_id()
				if did in dda:
					ld = dda[did]
				
				# check for acc
				aid = junc.get_acc_id()
				if aid in dda:
					la = dda[aid]
				
				# annotate the junction
				if len(la) > 0 and len(ld) > 0:
					# both ends are annotated - use an intersection
					junc.tid = list(set(ld[0] + la[0])) + ["*"]
					junc.gid = list(set(ld[1] + la[1]))
					junc.gname = list(set(ld[2] + la[2]))
					junc.strand = ld[3]
				else:
					if len(ld) > 0:
						# copy donor info
						junc.tid = list(ld[0]) + ["*"]
						junc.gid = list(ld[1])
						junc.gname = list(ld[2])
						junc.strand = ld[3]
					elif len(la) > 0:
						# copy acceptor info
						junc.tid = list(la[0]) + ["*"]
						junc.gid = list(la[1])
						junc.gname = list(la[2])
						junc.strand = la[3]
				
	#
	# next step is to produce a GTF for all of the junctions where each transcript
	# represents a single junction or a junction boundary (i.e. a theta quantification)
	#
	
	# open gtf file
	fout = open(juncdb_gtf, "w")
	sys.stderr.write("writing juncdb annotation {}\n".format(juncdb_gtf))
	for jid in djuncs.keys():
		junc = djuncs[jid]
		
		# per junction we'll make 4 GTF rows. two for the splice and two for the 
		# 5' and 3' retention
		
		grow = GtfRow(junc.rname, "psi", "exon", junc.start-args.l, junc.start-1, junc.strand)
		# set attributes
		
		grow.tid = junc.id
		# are we annotated?
		if len(junc.tid) > 0:
			# yes
			grow.attrs["oId"] = ",".join(junc.tid)
			grow.attrs["gene_name"] = ",".join(junc.gname)
			grow.gid = ",".join(junc.gid)
		else:
			# no - make something up
			ngid = "NOVG_{:08d}".format(ngidx)
			ngidx += 1
			grow.gid = ngid
			grow.attrs["gene_name"] = ngid
			grow.attrs["oId"] = ngid
		
		# ok, good to go
		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")
		
		# modify to make the other end
		grow.start = junc.end+1
		grow.end = junc.end+args.l
		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")
		
		# modify to make the theta rows
		grow.db = "theta"
		grow.start = junc.start-args.l
		grow.end = junc.start+args.l
		grow.tid = junc.get_donor_id()
		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")
		
		grow.start = junc.end-args.l
		grow.end = junc.end+args.l
		grow.tid = junc.get_acc_id()
		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")
	
	fout.close()
	sys.stderr.write("generating FASTA {}\n".format(juncdb_fa))
	cmd = "gffread -g {} -w {} {}".format(args.fa, juncdb_fa, juncdb_gtf)
	runcmd(cmd)

	# kmer the reference
	cmd = "kmercountexact.sh in={} k={} tuc=t out=temp_mers.fa overwrite=t".format(juncdb_fa, args.l)
	runcmd(cmd)
	
	# match back to the reference with seal
	cmd = "seal.sh in=temp_mers.fa ref={} k={} rcomp=t rename=t tuc=t hdist=0 mm=f trd=t out=temp_mers_matched.fa overwrite=t".format(juncdb_fa, args.l)
	runcmd(cmd)
	
	# now run that "other" python script to make the final index
	cmd = "{}/coding/python/juncdb-filter-matched-ref-mers.py {} {}".format(HOME, juncdb_gtf, "temp_mers_matched.fa")
	# open output file
	fout = open(juncdb_index, "w")
	# run the command 
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split(), stdout=fout)
	p1.wait()
	fout.close()
	
	# clean up
	os.unlink("temp_mers.fa")
	os.unlink("temp_mers_matched.fa")
	
	# compress 
	runcmd("gzip {}".format(juncdb_fa))
	runcmd("gzip {}".format(juncdb_index))

	return 0


#==============================================================================
# classes
#==============================================================================


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


# class for junction objects
class Junction(object):
	def __init__(self, rname="", start=0, end=0, id="", index=0, 
				strand="+", count=0, tid=None, gname=None, gid=None):
		
		self.rname = rname
		self.start = int(start)
		self.end = int(end)
		self.id = id
		self.index = index
		self.strand = strand
		self.count = float(count)
		
		if tid is None:
			self.tid = []
		else:
			self.tid = list(tid)
		
		if gname is None:
			self.gname = []
		else:
			self.gname = list(gname)	

		if gid is None:
			self.gid = []
		else:
			self.gid = list(gid)
	
	# --
	# compare this to another junction object
	def compare(self, jobj):
		res = False
		
		if self.rname==jobj.rname and self.start==jobj.start and self.end==jobj.end:
			res = True
		
		return res
	
	#-- 
	# merge a second Junction with the same position values into this one. this
	# just collapses in the annotation info
	def merge(self, jobj):
		# not sure what to do with index when we merge these
		self.index = -1
		
		# update tid, gid and gname
		tmp = set(self.tid)
		tmp.update(jobj.tid)
		self.tid = list(tmp)
		
		tmp = set(self.gid)
		tmp.update(jobj.gid)
		self.gid = list(tmp)
		
		tmp = set(self.gname)
		tmp.update(jobj.gname)
		self.gname = list(tmp)
		
		return 0
	
	def set_id_from_position(self):
		self.id = "{}:{}-{}".format(self.rname, self.start, self.end)
		return 0
	
	# return string of the donor id
	def get_donor_id(self):
		did = "{}:{}:5p".format(self.rname, self.start)
		return did
	
	# return string of the acceptor id
	def get_acc_id(self):
		aid = "{}:{}:3p".format(self.rname, self.end)
		return aid
			
			

#==============================================================================
# functions
#==============================================================================

# 
# parse_juncs
# parse a junctions file into a dict
def parse_juncs(f, skip_lines):
	# parse 
	
	# variables
	szl = ""
	aln = []
	junc = None
	ljuncs = []
	
	# open file
	fin = open(f, "r")
	
	if skip_lines > 0:
		for i in range(skip_lines):
			szl = fin.readline()
	
	for szl in fin:
		aln = szl.strip().split("\t")
		
		junc = Junction(rname=aln[0], start=aln[1], end=aln[2])
		junc.set_id_from_position()
		ljuncs.append(copy.deepcopy(junc))
	
	fin.close()
	
	return(ljuncs)

#--
# parse_gtf_to_junctions
# creates a list of Junction objects
# calling code should verify that the file exists and maybe throw this
# whole thing in a 'try'
def parse_gtf_to_junctions(f):
	
	ljuncs = []
	ltmp = []
	tid = ""
	last_tid = ""
	lrow = []
	aln = []
	szl = ""
	jidx = 0
	junc = None
	
	# open file
	fin = open(f, "r")
	
	# loop through gtf
	for szl in fin:
		aln = szl.strip().split("\t")
		attr = parse_gtf_attr(aln[8])
		tid = attr['transcript_id']
		
		if tid == last_tid:
			# make junction
			
			# keep index of junction
			jidx += 1
			# make a new Junction
			junc = Junction(rname=aln[0], start=int(lrow[4])+1, end=int(aln[3])-1, index=jidx, strand=aln[6])
			# make an id so that i don't have to do it later
			junc.set_id_from_position()
			# append in the id and the attributes
			junc.tid.append(attr['transcript_id'])
			junc.gid.append(attr['gene_id'])
			if "gene_name" in attr:
				junc.gname.append(attr['gene_name'])
				
			# add junction to the list - clearly there will be redundant junctions as we go
			ljuncs.append(copy.deepcopy(junc))
			
		else:
			jidx = 0
			
		lrow = list(aln)
		last_tid = tid 
	
	fin.close()
	
	return(ljuncs)


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
# junc_list_to_dict
# From a list of junctions to a dict keyed by the junction ids. matching junctions are 
# merged. annotation will propagate into junctions in this way.
def junc_list_to_dict(ljuncs):
	djuncs = {}
	n = len(ljuncs)
	
	# loop through
	for i in range(n):
		if ljuncs[i].id not in djuncs:
			djuncs[ljuncs[i].id] = copy.deepcopy(ljuncs[i])
		
		else:
			djuncs[ljuncs[i].id].merge(ljuncs[i])
	
	return djuncs

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


parser = argparse.ArgumentParser(description="Build junction database for mapping.")
parser.add_argument('stub', type=str, help="Name stub for the build")
parser.add_argument('fa', type=str, help="Genome FASTA reference")
parser.add_argument('juncs', type=str, nargs="+", help="Junctions file(s)")
parser.add_argument("-g", "--gtf", type=str, default=None, 
				help="GTF annotation to include known junctions and provide annotation")
parser.add_argument('-o', type=int, action="store", default=4, 
	help="Base pairs of overlap you'll want to require for valid hits [4]")
parser.add_argument('-l', type=int, action="store", default=26, 
	help="Read length [26]")
parser.add_argument('-s', type=int, default=0, 
				help="Skip this many lines from top of junctions files [0]")

args = parser.parse_args()

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

