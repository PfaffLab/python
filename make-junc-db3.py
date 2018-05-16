#!/usr/bin/env python
#==============================================================================
# make-junc-db3.py
#
# Shawn Driscoll
# 20160425
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Trying to make this simple....
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

	#
	# load junctions files
	#

	djuncs = {}
	dannot = {}
	dda = {}
	daa = {}

	for i in range(len(args.juncs)):
		sys.stderr.write("loading junctions from {}\n".format(args.juncs[i]))
		fin = open(args.juncs[i], "r")

		# skip lines?
		if args.s > 0:
			j = 0
			while j < args.s:
				sin = fin.readline()
				j += 1

		# loop through junctions in file
		for szl in fin:
			aln = szl.strip().split("\t")
			jid, lid, rid = pos_to_jid_set(aln[0:3])

			# add junction into dict
			if jid not in djuncs:

				djuncs[jid] = Junc(chrom=aln[0], start=int(aln[1]), end=int(aln[2]))
				djuncs[jid].update_id()

			# add left (donor) end into left side dict
			if lid not in dda:
				dda[lid] = JuncDA(chrom=aln[0], pos=int(aln[1]), parent=jid)
			else:
				dda[lid].update_parent(jid)

			# add right (acceptor) end into right side dict
			if rid not in daa:
				daa[rid] = JuncDA(chrom=aln[0], pos=int(aln[2]), parent=jid)
			else:
				daa[rid].update_parent(jid)

		# done close file
		fin.close()

	# finished loading junctions - is there an annotation?
	if args.gtf is not None:
		# load annotation
		sys.stderr.write("loading annotation from {}\n".format(args.gtf))
		fin = open(args.gtf, "r")

		for szl in fin:

			grow = GtfRow()
			grow.parse(szl.strip())

			if grow.type=="exon":

				if grow.tid not in dannot:
					dannot[grow.tid] = []

				dannot[grow.tid].append(grow)

		fin.close()

		# loaded now we have to find all of the junctions
		ltid = dannot.keys()
		sys.stderr.write("loaded {} transcripts. parsing out junctions and annotation\n".format(len(ltid)))

		for tid in ltid:
			if len(dannot[tid]) > 1:
				# this is a multi-exon feature. make sure the features are sorted
				lstarts = []
				lends = []
				for i in range(len(dannot[tid])):
					lstarts.append(dannot[tid][i].start)
					lends.append(dannot[tid][i].end)

				o = list(np.argsort(lstarts))

				for i in range(1, len(o)):
					# junction is from this position's start to the last position's end
					growA = dannot[tid][i-1]
					growC = dannot[tid][i]

					jtmp = Junc(chrom=growA.rname, start=growA.end+1, end=growC.start-1)
					jtmp.update_id()
					jtmp.strand = growA.strand

					if jtmp.jid not in djuncs:
						# this junction is not in the database, add it
						djuncs[jtmp.jid] = jtmp

					if jtmp.lid not in dda:
						dda[jtmp.lid] = JuncDA(chrom=jtmp.chrom, pos=jtmp.start, parent=jtmp.jid)
					else:
						dda[jtmp.lid].update_parent(jtmp.jid)

					if jtmp.rid not in daa:
						daa[jtmp.rid] = JuncDA(chrom=jtmp.chrom, pos=jtmp.end, parent=jtmp.jid)
					else:
						daa[jtmp.rid].update_parent(jtmp.jid)

					# the junction was either there already or we just added it so now we can 
					# update its annotation

					djuncs[jtmp.jid].annotated = True
					djuncs[jtmp.jid].update_annotation(growA.tid, growA.gid, growA.attrs['gene_name'])

					dda[jtmp.lid].annotated = True
					dda[jtmp.lid].update_annotation(growA.tid, growA.gid, growA.attrs['gene_name'])

					daa[jtmp.rid].annotated = True
					daa[jtmp.rid].update_annotation(growA.tid, growA.gid, growA.attrs['gene_name'])

		# finished with that, now we have to propagate annotation into all of the junctions via the 
		# donor/acceptor lists

		for jid in djuncs.keys():
			
			if not djuncs[jid].annotated:
				# this one isn't annotated - check the donor and acceptor
				dlid = dda[djuncs[jid].lid]
				drid = daa[djuncs[jid].rid]

				# the junction isn't annotated so therefore there is not a transcript id for it
				(djuncs[jid].transcript_id).update("*")

				# figure out gene id, gene name and strand
				if dlid.annotated and drid.annotated:
					
					djuncs[jid].annotated = True

					# find strand
					for pjid in list(dlid.parent):
						if djuncs[pjid].annotated:
							djuncs[jid].strand = djuncs[pjid].strand
							break

					sgid = (dlid.gene_id).intersection(drid.gene_id)
					if len(sgid)==0:
						# no gene id shares the donor and acceptor although they 
						# are both annotated. i think we should drop this one
						djuncs[jid].keep = False
					else:
						# use intersection. sometimes there are more than 1 gene id because of the annotation
						# which sometimes has multiple gene ids in a single gene locus
						(djuncs[jid].gene_id).update(list(sgid))

					sgname = (dlid.gene_name).intersection(drid.gene_name)
					if len(sgname)==0:
						(djuncs[jid].gene_name).update("*")
					else:
						(djuncs[jid].gene_name).update(list(sgname))

				else:

					# left side is annotated
					if dlid.annotated:
						djuncs[jid].annotated = True

						# find strand
						for pjid in list(dlid.parent):
							if djuncs[pjid].annotated:
								djuncs[jid].strand = djuncs[pjid].strand
								break

						# copy gene id and gene name
						(djuncs[jid].gene_id).update(list(dlid.gene_id))
						(djuncs[jid].gene_name).update(list(dlid.gene_name))

					elif drid.annotated:
						djuncs[jid].annotated = True

						# find strand
						for pjid in list(drid.parent):
							if djuncs[pjid].annotated:
								djuncs[jid].strand = djuncs[pjid].strand
								break

						# copy gene id and gene name
						(djuncs[jid].gene_id).update(list(drid.gene_id))
						(djuncs[jid].gene_name).update(list(drid.gene_name))


				
	#
	# next step is to produce a GTF for all of the junctions where each transcript
	# represents a single junction or a junction boundary (i.e. a theta quantification)
	#
	
	# open gtf file
	fout = open(juncdb_gtf, "w")
	sys.stderr.write("writing juncdb annotation {}\n".format(juncdb_gtf))
	for jid in djuncs.keys():
		junc = djuncs[jid]
		if not junc.keep:
			continue
		
		# per junction we'll make 2 GTF rows
		
		grow = GtfRow(junc.chrom, "psi", "exon", junc.start-args.l, junc.start-1, junc.strand)
		# set attributes
		
		grow.tid = junc.jid
		# are we annotated?
		if junc.annotated:
			# yes
			grow.attrs["oId"] = ",".join(junc.transcript_id)
			grow.attrs["gene_name"] = ",".join(junc.gene_name)
			grow.gid = ",".join(junc.gene_id)
		else:
			# no - make something up
			ngid = "NOVG_{:08d}".format(ngidx)
			ngidx += 1
			grow.gid = ngid
			grow.attrs["gene_name"] = ngid
			grow.attrs["oId"] = ngid
			# update the donor and acceptor for this junction so these names are 
			# propagated into the gtf
			dda[junc.lid].update_annotation(ngid, ngid, ngid)
			daa[junc.rid].update_annotation(ngid, ngid, ngid)
		
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

	# make theta features
	for lid in dda.keys():
		dlid = dda[lid]

		grow = GtfRow(rname=dlid.chrom, db="theta", type="exon", start=dlid.pos-args.l, end=dlid.pos+args.l, strand="+")
		grow.tid = "{}:{}".format(dlid.chrom, dlid.pos)

		# if annotated then get the strand otherwise assume it's positive
		if dlid.annotated:
			for pjid in list(dlid.parent):
				if djuncs[pjid].annotated:
					grow.strand = djuncs[pjid].strand
					break

		# update attributes
		grow.attrs["oId"] = ",".join(dlid.transcript_id)
		grow.attrs["gene_name"] = ",".join(dlid.gene_name)
		grow.gid = ",".join(dlid.gene_id)

		if grow.strand == "+":
			grow.tid += ":5p"
		else:
			grow.tid += ":3p"

		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")

	# make theta features
	for rid in daa.keys():
		drid = daa[rid]

		grow = GtfRow(rname=drid.chrom, db="theta", type="exon", start=drid.pos-args.l, end=drid.pos+args.l, strand="+")
		grow.tid = "{}:{}".format(drid.chrom, drid.pos)

		# if annotated then get the strand otherwise assume it's positive
		if drid.annotated:
			for pjid in list(drid.parent):
				if djuncs[pjid].annotated:
					grow.strand = djuncs[pjid].strand
					break

		# update attributes
		grow.attrs["oId"] = ",".join(drid.transcript_id)
		grow.attrs["gene_name"] = ",".join(drid.gene_name)
		grow.gid = ",".join(drid.gene_id)

		if grow.strand == "+":
			grow.tid += ":3p"
		else:
			grow.tid += ":5p"

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

class Junc(object):
	def __init__(self, chrom="", start=0, end=0):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.strand = "+"
		self.annotated = False
		self.jid = ""
		self.lid = ""
		self.rid = ""
		self.gene_id = set()
		self.transcript_id = set()
		self.gene_name = set()
		self.keep = True
		return None

	def update_id(self):
		self.jid = "{}:{}-{}".format(self.chrom, self.start, self.end)
		self.lid = "{}:{}".format(self.chrom, self.start)
		self.rid = "{}:{}".format(self.chrom, self.end)
		return 0

	def update_annotation(self, tid, gid, gname):
		
		if type(tid) is not list:
			tid = [tid]
		if type(gid) is not list:
			gid = [gid]
		if type(gname) is not list:
			gname = [gname]

		(self.transcript_id).update(tid)
		(self.gene_id).update(gid)
		(self.gene_name).update(gname)
		return 0

class JuncDA(object):
	def __init__(self, chrom="", pos=0, parent=None):
		self.chrom = chrom
		self.pos = pos
		
		if parent is not None:
			self.parent = set([parent])
		else:
			self.parent = set()

		self.annotated = False
		self.gene_id = set()
		self.transcript_id = set()
		self.gene_name = set()
		return None

	def update_parent(self, pid):
		if type(pid) is not list:
			pid = [pid]
		
		(self.parent).update(pid)
		return 0

	def update_annotation(self, tid, gid, gname):
		
		if type(tid) is not list:
			tid = [tid]
		if type(gid) is not list:
			gid = [gid]
		if type(gname) is not list:
			gname = [gname]

		(self.transcript_id).update(tid)
		(self.gene_id).update(gid)
		(self.gene_name).update(gname)
		return 0

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

def pos_to_jid_set(pos):

	jid = "{}:{}-{}".format(pos[0], pos[1], pos[2])
	lid = "{}:{}".format(pos[0], pos[1])
	rid = "{}:{}".format(pos[0], pos[2])
	return jid, lid, rid



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

