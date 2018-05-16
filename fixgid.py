#!/usr/bin/python
#
# fixgid.py
#
# This processes a gtf file with both gene_id, gene_name and locus attributes
# and makes a new set of gene ids from all unique combinaions of gene_name and 
# locus that exist in the GTF. Then a second pass combines all genes with the 
# exact same name into gene id clusters
#

import sys, re
from hashlib import md5
import argparse
import igraph as ig


parser = argparse.ArgumentParser(description="Replace gene_id field in GTF by generating a unique ID per combination of gene name and locus id.")
parser.add_argument("gtf", type=str, help="GTF file to modify")
parser.add_argument('-s', action="store", dest="stub", type=str, default="GID", 
	help="Stub for gene_id (default: GID)")
parser.add_argument('-m', '--mode', action="store", type=str, default="log", 
	help="Mode for grouping features into 'gene_id'. Default is to use 'locus' or 'gene_name' attributes (log). You may also use 'locus' and 'gene_name' (lag), just 'locus' (l) or just 'gene_name' (g)")

args = parser.parse_args()

# check modes
modes = set(['log', 'lag', 'l', 'g'])
if args.mode not in modes:
	sys.stderr.write("Error: unknown mode provided with -m. Mode must be 'log', 'lag', 'l', or 'g'\n")
	sys.exit(1)

# --
# parse the attributes field of a GTF row into a hash
def parse_gtf_attr(sz):
	# split string on ;
	s1 = sz.split(";")
	attr = {}

	for i in range(len(s1)):

		tmp = re.sub("^\s+", "", s1[i])
		tmp = re.sub("\s+$", "", tmp)

		if len(s1[i]) > 0:
			# split by space
			s2 = tmp.split(" ")
			tmp = re.sub("^\"", "", s2[1])
			tmp = re.sub("\"$", "", tmp)

			attr[s2[0]] = tmp

	return attr

def message(sz):
	sys.stderr.write("[fixgid] " + sz + "\n")

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

	#
	# this function updates the gene id of a row object
	def update_gid(self, gid_new):
		self.gid = gid_new
		self.attrs["gene_id"] = gid_new

		return 0


#==============================================================================
#==============================================================================
#
# BEGIN MAIN SCRIPT
#
#==============================================================================
#==============================================================================


# variables

dannot = {}
# list to maintain order of transcripts in the file
tid_order = []

# keep three hashes to connect locus and gene name to transcript id
l2tid = {}
gn2tid = {}
lgn2tid = {}
tid2idx = {}
tidx = 0

# open the GTF , parse stuff and build the hash of name/locus combinations

message("loading annotation")
fin = open(args.gtf, "r")

for szl in fin:

	# make GtfRow object
	grow = GtfRow()
	# parse line from file
	grow.parse(szl)

	# insert into the dannot hash
	if grow.tid not in dannot:
		tid_order.append(grow.tid)
		dannot[grow.tid] = []
		tid2idx[grow.tid] = tidx
		# transcript counter
		tidx += 1

	dannot[grow.tid].append(grow)

	# locus
	loc = grow.attrs['locus'].lower()
	if loc not in l2tid:
		l2tid[loc] = set()

	l2tid[loc].update([grow.tid])

	# gene name
	gn = grow.attrs['gene_name'].lower()
	if gn not in gn2tid:
		gn2tid[gn] = set()

	gn2tid[gn].update([grow.tid])

	lgn = "{}.{}".format(loc, gn)
	if lgn not in lgn2tid:
		lgn2tid[lgn] = set()

	lgn2tid[lgn].update([grow.tid])

fin.close()

#
# build edge table
message("building edge table")
ledges = []
if args.mode=="log" or args.mode=="l":
	# loop through the locus to tid hash and build edges
	for lid in l2tid.keys():
		tmp = list(l2tid[lid])
		
		for i in range(len(tmp)):
			for j in range(len(tmp)):
				ledges.append((tid2idx[tmp[i]], tid2idx[tmp[j]]))

if args.mode=="log" or args.mode=="g":
	# loop through 'gene name' to 'tid' hash
	for gn in gn2tid.keys():
		tmp = list(gn2tid[gn])

		for i in range(len(tmp)):
			for j in range(len(tmp)):
				ledges.append((tid2idx[tmp[i]], tid2idx[tmp[j]]))

if args.mode=="lag":
	# loop through the locus/gn hash
	for lgn in lgn2tid.keys():
		tmp = list(lgn2tid[lgn])
		for i in range(len(tmp)):
			for j in range(len(tmp)):
				ledges.append((tid2idx[tmp[i]], tid2idx[tmp[j]]))

# build graph and then cluster
message("building and clustering graph")
g = ig.Graph()
g.add_vertices(tidx)
g.add_edges(ledges)

tid_groups = g.clusters().membership

message("generating gene_id from clustering and writing to STDOUT")
for i in range(len(tid_groups)):
	gid_new = "{}_{:08d}".format(args.stub, tid_groups[i])
	tid = tid_order[i]

	for grow in dannot[tid]:
		grow.update_gid(gid_new)
		lout = grow.tolist()
		print "\t".join(map(str, lout))




