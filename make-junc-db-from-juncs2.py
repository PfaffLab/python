#!/usr/bin/env python
#==============================================================================
# make-junc-db-from-juncs2.py
#
# Shawn Driscoll
# 20160425
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This one is a very good version but cannot account for possible issues in the 
# annotation (provided with -g). I have found in my annotation sometimes 
# anti-sense transcripts that share many of their introns and splice sites which
# confuses the strand information when building a juncdb GTF. since this is a
# novel assembly anyways I'm making it so, for the output, everything is positive
# strand and any quantification against it should be unstranded.
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
	dgid_strand = {}

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

			# check junction length. you have to keep an eye on those observed junctions from alignments
			# because they get a little out of control sometimes.
			if int(aln[2])-int(aln[1]) > args.max_junc_len:
				continue

			# add junction into dict
			if jid not in djuncs:

				djuncs[jid] = Junc(chrom=aln[0], start=int(aln[1]), end=int(aln[2]))
				djuncs[jid].update_id()

# REMOVING STRAND STUFF FOR NOW
#				djuncs[jid].strand = "+"
			
			# add hit
			djuncs[jid].add_hit(aln[args.count_col])

		# done close file
		fin.close()

	sys.stderr.write("loaded {} junctions\n".format(len(djuncs.keys())))

	# filter the junctions out prior to moving forward...
	if args.min_count > 0:
		sys.stderr.write("filtering junctions to min count depth of {} in min {} sample(s)...\n".format(args.min_count, args.min_samples))

		# make a list of the jid's to drop
		jdrops = []
		for jid in djuncs.keys():

			n = djuncs[jid].num_counts_at_depth(args.min_count)
			if n < args.min_samples:
				jdrops.append(jid)

		if len(jdrops) > 0:
			sys.stderr.write("dropping {} junctions that did not meet minimum depth/samples...\n".format(len(jdrops)))
			for jid in jdrops:
				del(djuncs[jid])

			sys.stderr.write("retained {} junctions\n".format(len(djuncs.keys())))


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

# REMOVING STRAND STUFF FOR NOW
#				# track strand for annotated gene ids
#				if grow.gid not in dgid_strand:
#					dgid_strand[grow.gid] = grow.strand

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

				# get sort order of the exons
				o = list(np.argsort(lstarts))

				# loop through sorted exons and make junctions
				for i in range(1, len(o)):
					# junction is from this position's start to the last position's end
					growA = dannot[tid][i-1]
					growC = dannot[tid][i]

					jtmp = Junc(chrom=growA.rname, start=growA.end+1, end=growC.start-1)
					jtmp.update_id()

					if jtmp.jid not in djuncs:
						# this junction is not in the database, add it
						djuncs[jtmp.jid] = jtmp

					# the junction was either there already or we just added it so now we can 
					# update its annotation

					djuncs[jtmp.jid].annotated = True
					djuncs[jtmp.jid].update_annotation(growA.tid, growA.gid, growA.attrs['gene_name'])
					# updatae strand
# REMOVING STRAND STUFF FOR NOW
#					djuncs[jtmp.jid].strand = growA.strand

		# print update to total junctions
		sys.stderr.write("annotation junctions + observed: {}\n".format(len(djuncs.keys())))

	# --
	#
	# build the donor and acceptor dicts
	#
	# --

	dda = {}
	daa = {}
	sys.stderr.write("building donor/acceptor tables\n")
	for jid in djuncs.keys():

		junc = djuncs[jid]

		if junc.lid not in dda:
			# add it
			dda[junc.lid] = JuncDA(chrom=junc.chrom, pos=junc.start, parent=junc.jid)
			dda[junc.lid].update_id()
		else:
			dda[junc.lid].update_parent(junc.jid)

		if junc.rid not in daa:
			# add it
			daa[junc.rid] = JuncDA(chrom=junc.chrom, pos=junc.end, parent=junc.jid)
			daa[junc.rid].update_id()
		else:
			daa[junc.rid].update_parent(junc.jid)

		if junc.annotated:
			# update annotation of the donor/acceptor
			dda[junc.lid].annotated = True
			dda[junc.lid].update_annotation(list(junc.transcript_id), list(junc.gene_id), list(junc.gene_name))
# REMOVING STRAND STUFF FOR NOW
#			dda[junc.lid].strand = junc.strand

			daa[junc.rid].annotated = True
			daa[junc.rid].update_annotation(list(junc.transcript_id), list(junc.gene_id), list(junc.gene_name))
# REMOVING STRAND STUFF FOR NOW
#			daa[junc.rid].strand = junc.strand

	
	# --
	#
	# if we have a GTF then at this point all of the junctions that were in the GTF are 
	# now anntated and all of the donors/acceptors that were a part of those junctions
	# are also annotated. some donor/acceptors are NOT annotated and some junctions are 
	# not annotated.  the unannotated junctions are of two types: using annotated donor
	# and/or acceptor or not. the 'not' ones need novel ids which then need to be
	# propagated to the donor and acceptor.
	#
	# --
				
	#
	# next step is to produce a GTF for all of the junctions where each transcript
	# represents a single junction or a junction boundary (i.e. a theta quantification)
	#
	
	# open gtf output file

	fout = open(juncdb_gtf, "w")
	sys.stderr.write("writing juncdb annotation {}\n".format(juncdb_gtf))
	ngidx = 0

	# use these dicts to track novel gene names created for novel junctions as assigned to 
	# donors and acceptors
	dda_novG = {}
	daa_novG = {}

	for jid in sorted(djuncs.keys()):
		# make a references to the junction
		junc = djuncs[jid]


		if args.gtf is not None and not junc.annotated:
			# check the donor and acceptor for annotation
			dlid = dda[junc.lid]
			drid = daa[junc.rid]

			# no matter what happens below this is a novel junction
			(junc.transcript_id).update(["*"])

			if dlid.annotated and drid.annotated:
				lgid = (dlid.gene_id).intersection(drid.gene_id)
				lgname = (dlid.gene_name).intersection(drid.gene_name)

				if len(lgid) > 0:
					# good, update it
					(junc.gene_id).update(list(lgid))
				else:
					# no intersection between the two annotated features....weird
					sys.stderr.write("warning: junction between different annotated genes {}\n".format(junc.jid))
					junc.keep = False

				(junc.gene_name).update(list(lgname))
# REMOVING STRAND STUFF FOR NOW
#				junc.strand = dlid.strand
				junc.annotated = True

			else:

				if dlid.annotated:
					# use annotation from left side
					junc.strand = dlid.strand
					(junc.gene_id).update(list(dlid.gene_id))
					(junc.gene_name).update(list(dlid.gene_name))
					junc.annotated = True

					# this means right side is not annotated so we can pass some info over
					drid.update_annotation("*", list(dlid.gene_id), list(dlid.gene_name))
# REMOVING STRAND STUFF FOR NOW
#					drid.strand = dlid.strand
					drid.annotated = True

				elif drid.annotated:
					# use annotation from right
					junc.strand = drid.strand
					(junc.gene_id).update(list(drid.gene_id))
					(junc.gene_name).update(list(drid.gene_name))
					junc.annotated = True

					# this means left side is not annotated so we can pass some info over
					dlid.update_annotation("*", list(drid.gene_id), list(drid.gene_name))
# REMOVING STRAND STUFF FOR NOW
#					dlid.strand = drid.strand
					dlid.annotated = True

				else:
					# neither end is annotated so we need a new id for this guy

					if dlid.id in dda_novG and drid.id in daa_novG:
						# both already have a novel id
						ngid = [dda_novG[dlid.id], daa_novG[drid.id]]
						if ngid[0] != ngid[1]:
							sys.stderr.write("warning: different novel gene ids at novel junction {} and {}\n".format(dlid.id, drid.id))
					else:
						# check if either of the two have a novel id already
						if dlid.id in dda_novG:
							# get id from left side
							ngid = dda_novG[dlid.id]
							# set it for the right side
							daa_novG[drid.id] = ngid
						elif drid.id in daa_novG:
							# get id from right side
							ngid = daa_novG[drid.id]
							# set it for the left side
							dda_novG[dlid.id] = ngid
						else:
							# make a new id						
							ngid = "NOVG_{:08d}".format(ngidx)
							ngidx += 1
							# set for both
							dda_novG[dlid.id] = ngid
							daa_novG[drid.id] = ngid

					# update all of the objects with ngid which should be defined
					if type(ngid) is not list:
						ngid = [ngid]

					# update everything
					dlid.update_annotation("*", ngid, ngid)
					drid.update_annotation("*", ngid, ngid)
					(junc.gene_id).update(ngid)
					(junc.gene_name).update(ngid)


		if not junc.keep:
			# next...
			continue

# REMOVING STRAND STUFF FOR NOW
#		if junc.annotated:
#			# check the strand just in case
#			jstrands = set()
#			for gid in list(junc.gene_id):
#				if gid in dgid_strand:
#					jstrands.update([dgid_strand[gid]])
#			jstrands = list(jstrands)
#			if len(jstrands)==0:
#				print junc.jid, junc.gene_id, len(junc.gene_id)
#				return 1
#			elif len(jstrands) > 1:
#				# what is going on here?
#				sys.stderr.write("warning: dropping junction {} with two strands @ gene id ({})\n".format(jid, ",".join(junc.gene_id)))
#				junc.keep = False
#			elif junc.strand != jstrands[0]:
#				sys.stderr.write("warning: updating incorrect strand at {}\n".format(jid))
#				junc.strand = jstrands[0]


		
		# per junction we'll make 2 GTF rows
		
		grow = GtfRow(junc.chrom, "psi", "exon", junc.start-args.l, junc.start-1, junc.strand)
		# set attributes
		
		grow.tid = junc.jid
		grow.attrs["oId"] = ",".join(junc.transcript_id)
		grow.attrs["gene_name"] = ",".join(junc.gene_name)
		grow.gid = ",".join(junc.gene_id)
		
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
	for lid in sorted(dda.keys()):
		dlid = dda[lid]

# REMOVING STRAND STUFF FOR NOW
#		if args.gtf is not None:
#			# check the strand
#			jstrands = set()
#			for gid in list(dlid.gene_id):
#				if gid in dgid_strand:
#					jstrands.update([dgid_strand[gid]])
#			jstrands = list(jstrands)
#			if len(jstrands)==0:
#				sys.stderr.write("error: no strand info on an annotated donor {}\n".format(dlid.id))
#				return 1
#			elif len(jstrands) > 1:
#				# what is going on here?
#				sys.stderr.write("warning: dropping donor {} with two strands @ gene id ({})\n".format(dlid.id, ",".join(dlid.gene_id)))
#				dlid.keep = False
#			elif dlid.strand != jstrands[0]:
#				sys.stderr.write("warning: updating incorrect strand at {}\n".format(dlid.id))
#				dlid.strand = jstrands[0]

		if not dlid.keep:
			continue

		grow = GtfRow(rname=dlid.chrom, db="theta", type="exon", start=dlid.pos-args.l, end=dlid.pos+args.l, strand=dlid.strand)
		grow.tid = "{}:{}".format(dlid.chrom, dlid.pos)

		# update attributes
		grow.attrs["oId"] = ",".join(dlid.transcript_id)
		grow.attrs["gene_name"] = ",".join(dlid.gene_name)
		grow.gid = ",".join(dlid.gene_id)

# REMOVING STRAND STUFF FOR NOW
#		if grow.strand == "+":
#			grow.tid += ":5p"
#		else:
#			grow.tid += ":3p"

		grow.tid += ":5p"

		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")

	# make theta features
	for rid in sorted(daa.keys()):
		drid = daa[rid]

# REMOVING STRAND STUFF FOR NOW
#		if args.gtf is not None:
#			# check the strand
#			jstrands = set()
#			for gid in list(drid.gene_id):
#				if gid in dgid_strand:
#					jstrands.update([dgid_strand[gid]])
#			jstrands = list(jstrands)
#			if len(jstrands)==0:
#				sys.stderr.write("error: no strand info on an annotated acceptor {}\n".format(drid.id))
#				return 1
#			elif len(jstrands) > 1:
#				# what is going on here?
#				sys.stderr.write("warning: dropping acceptor {} with two strands @ gene id ({})\n".format(drid.id, ",".join(drid.gene_id)))
#				drid.keep = False
#			elif drid.strand != jstrands[0]:
#				sys.stderr.write("warning: updating incorrect strand at {}\n".format(drid.id))
#				drid.strand = jstrands[0]

		if not drid.keep:
			continue

		grow = GtfRow(rname=drid.chrom, db="theta", type="exon", start=drid.pos-args.l, end=drid.pos+args.l, strand=drid.strand)
		grow.tid = "{}:{}".format(drid.chrom, drid.pos)

		# update attributes
		grow.attrs["oId"] = ",".join(drid.transcript_id)
		grow.attrs["gene_name"] = ",".join(drid.gene_name)
		grow.gid = ",".join(drid.gene_id)

# REMOVING STRAND STUFF FOR NOW
#		if grow.strand == "+":
#			grow.tid += ":3p"
#		else:
#			grow.tid += ":5p"
		
		grow.tid += ":3p"

		lout = grow.tolist()
		fout.write("\t".join(map(str, lout)))
		fout.write("\n")
		
	# close output file
	fout.close()

	# --
	#
	# generate FASTA from the GTF and user specified genome sequence
	#
	# --

	sys.stderr.write("generating FASTA {}\n".format(juncdb_fa))
	cmd = "gffread -g {} -w {} {}".format(args.fa, juncdb_fa, juncdb_gtf)
	runcmd(cmd)

	# --
	#
	# kmer the FASTA
	#
	# --

	cmd = "kmercountexact.sh in={} k={} tuc=t out=temp_mers.fa overwrite=t".format(juncdb_fa, args.l)
	runcmd(cmd)
	
	# --
	#
	# match back to the reference with seal
	#
	# --

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
		self.counts = []
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

	def add_hit(self, v):
		(self.counts).append(float(v))
		return 0

	def num_counts_at_depth(self, d):
		n = len(self.counts)
		nc = 0

		for i in range(n):
			if self.counts[i] >= d:
				nc += 1

		return nc

class JuncDA(object):
	def __init__(self, chrom="", pos=0, parent=None):
		self.chrom = chrom
		self.pos = pos
		self.id = ""

		if parent is not None:
			self.parent = set([parent])
		else:
			self.parent = set()

		self.annotated = False
		self.gene_id = set()
		self.transcript_id = set()
		self.gene_name = set()
		self.strand = "+"
		self.keep = True

		return None

	def update_id(self):
		self.id = "{}:{}".format(self.chrom, self.pos)
		return 0

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

		return None
	
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
parser.add_argument('-j', '--max-junc-len', type=int, default=100000, 
	help="Maximum length for observed junctions (does not apply to those loaded from annotation [100000]")
parser.add_argument('-m', '--min-count', type=int, default=8, 
	help="Minimum count level for a junction to be considered from the junctions files. [8]")
parser.add_argument('-p', '--min-samples', type=int, default=1, 
	help="Minimum number of samples at minimum count depth (set with -m) [1]")
parser.add_argument('-c', '--count-col', type=int, default=3, 
	help="Zero-based column index for counts of junctions in the junctions files. Default assumes the psi/theta format [3]")
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

