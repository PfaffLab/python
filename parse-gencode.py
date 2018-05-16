#!/usr/bin/python
#
# parse-gencode.py
# 20170309
# Shawn Driscoll
#
#

import sys, re
import argparse


parser = argparse.ArgumentParser(description="Parse gencode annotation into a TSV")
parser.add_argument("gtf", type=str, help="Gencode GTF")
parser.add_argument("--tsv", action="store_const", const=True, default=False, 
	help="Export as TSV with a header. Otherwise GTF is returned")
parser.add_argument("-f", type=str, choices=["transcript", "gene", "exon"], default="transcript", 
	help="Select the type of feature to parse. May be 'transcript', 'gene', or 'exon' [transcript]")
parser.add_argument('-t', type=str, default="all", 
	help="Specify the gene type to extract (such as 'protein_coding'). You may specify more than one as a comma separated list with no spaces [all]")

args = parser.parse_args()

def load_gtf(fname, type, gtypes):
	dannot = {}
	dgid2tid = {}
	dtid2gid = {}
	dtid2gn = {}
	dgn2tid = {}
	dattrs = set()

	grow = None
	nrows = 0
	nparsed = 0

	fin = open(fname, "r")

	for szl in fin:

		if re.search("^#", szl):
			continue

		nrows += 1

		if (nrows % 100000) == 0:
			message("read: {}; parsed: {}".format(nrows, nparsed))

		aln = szl.split("\t")
		if aln[2] != type:
			continue

		grow = GtfRow()
		grow.parse(szl)

		if "all" not in gtypes:
			if "gene_type" in grow.attrs:
				if grow.attrs['gene_type'] not in gtypes:
					continue
			elif "transcript_type" in grow.attrs:
				if grow.attrs['transcript_id'] not in gtypes:
					continue

		nparsed += 1

		gname = grow.attrs['gene_name']
		tid = grow.tid

		if tid not in dannot:
			dannot[tid] = []
			dtid2gid[tid] = ""
			dtid2gn[tid] = ""

		dannot[tid].append(grow)
		dtid2gid[tid] = grow.gid
		dtid2gn[tid] = gname

		if gname not in dgn2tid:
			dgn2tid[gname] = set()

		dgn2tid[gname].update([tid])

		if grow.gid not in dgid2tid:
			dgid2tid[grow.gid] = set()

		dgid2tid[grow.gid].update([tid])

		dattrs.update(grow.attrs.keys())


	fin.close()

	return dannot, dgid2tid, dtid2gid, dtid2gn, dgn2tid, dattrs

def message(sz):
	sys.stderr.write("[parse-gencode] " + sz + "\n")

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

	def __str_(self):
		ll = self.tolist()
		return "\t".join(map(str, ll))
		
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
		fsplit = aln[8].split(";")
		for i in range(len(fsplit)):
			# remove leading and trailing whitespace
			fsplit[i] = re.sub("^\s+", "", fsplit[i])
			fsplit[i] = re.sub("\s+$", "", fsplit[i])
			
			# get key
			r = re.search("^([^\s]+)", fsplit[i])
			if r:
				kid = r.group(1)
				self.attrs[kid] = ""

				# get value
				if re.search("\"", fsplit[i]):
					# quoted value exists
					r = re.search("\"([^\"]+)\"", fsplit[i])
					if r:
						self.attrs[kid] = r.group(1)
				else:
					# unquoted value
					tmp = fsplit[i].split(" ")
					self.attrs[kid] = tmp[-1]

		if 'transcript_id' in self.attrs:
			self.tid = self.attrs['transcript_id']
		if 'gene_id' in self.attrs:
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
annot = {}
tid2gn = {}
gid2tid = {}
tid2gid = {}
gn2tid = {}
attrs = {}

# list to maintain order of transcripts in the file
tid_order = []

# keep three hashes to connect locus and gene name to transcript id
tidx = 0

# load gtf
annot, gid2tid, tid2gid, tid2gn, gn2tid, attrs = load_gtf(args.gtf, args.f, args.t)

primary_attrs = ["transcript_id", "gene_id", "gene_name"]
attrs = attrs - set(primary_attrs)
attrs = sorted(list(attrs))
attrs = primary_attrs + attrs

if args.tsv:
	header = ["ref", "db", "type", "start", "end", "foo", "strand", "bar"] + attrs
	print "\t".join(header)

for tid in sorted(annot.keys()):
	for grow in annot[tid]:
		lout = grow.tolist()

		if args.tsv:
			lout = lout[0:(len(lout)-1)]
			
			for k in attrs:
				if k in grow.attrs:
					lout.append(grow.attrs[k])
				else:
					lout.append("u")		

		print "\t".join(map(str, lout))

sys.exit(0)




