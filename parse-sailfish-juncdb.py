#!/usr/bin/env python
#==============================================================================
# parse-sailfish-juncdb.py
#
# Shawn Driscoll
# 20160825
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This file is meant to run on the output of Sailfish when quantifying 
# psi/theta junction database hits from raw reads. 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

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
	dgtf = None
	djdb = {}
	dtrans = {}
	dquant = {}

	# load the GTF annotation
	sys.stderr.write("Loading {}\n".format(args.jdb_gtf))
	rres = parse_gtf(args.jdb_gtf, djdb, dtrans)

	sys.stderr.write("Parsing {}\n".format(args.sfish_quant))

	# now parse the quantificaiton into this thing
	rres = parse_sfish(args.sfish_quant, djdb, dtrans)

	# run through djdb and produce output
	print jdb_header()
	for jid in sorted(djdb.keys()):
		# get total donor count
		did = djdb[jid]['did']
		if did in dtrans:
			djdb[jid]['donor_use'] = dtrans[did]['use_count']
			djdb[jid]['donor_ovl'] = dtrans[did]['skip_count']
		aid = djdb[jid]['aid']
		if aid in dtrans:
			djdb[jid]['acc_use'] = dtrans[aid]['use_count']
			djdb[jid]['acc_ovl'] = dtrans[aid]['skip_count']

		dout = djdb[jid]
		rres = update_metrics(dout)
		print jdb_tostring(dout)

	return 0



#==============================================================================
#
# function definitions
#
#==============================================================================

def update_metrics(d):

	if d['donor_use'] > 0 or d['acc_use'] > 0:
		# this means we can have psi and theta values. if neither the donor or
		# acceptor is used at all then this will keep the NA values for the
		# splice metrics

		if d['donor_use'] > 0:
			d['psi5'] = d['count']*1.0/(d['donor_use'])
			d['theta5'] = d['donor_use']*1.0/(d['donor_use']+d['donor_ovl'])

		if d['acc_use'] > 0:
			d['psi3'] = d['count']*1.0/(d['acc_use'])
			d['theta3'] = d['acc_use']*1.0/(d['acc_use']+d['acc_ovl'])

	return 0

def jdb_header():
	lout = ['chrom', 'start', 'end', 'count', 'id', 'strand', 'donor_ovl', 'donor_use', 
		'acc_ovl', 'acc_use', 'psi5', 'psi3', 'theta5', 'theta3', 'gene_name', 'gene_id', 'transcript_id']
	return "\t".join(lout)

def jdb_tostring(d):

	lout = [d['chrom'], d['start'], d['end'], d['count'], d['jid'], 
		d['strand'], d['donor_ovl'], d['donor_use'], d['acc_ovl'], d['acc_use'], 
		d['psi5'], d['psi3'], d['theta5'], d['theta3'], d['gene_name'], d['gene_id'], 
		d['transcript_id']]

	sz = "\t".join(map(str, lout))
	return sz


# 
# is_da
# this function returns true of the passed string is a donor/acceptor name
def is_da(sz):

	if re.search("[53]p$", sz):
		return True
	
	return False

#
# parse Sailfish output.
# f - input file name
# djuncs - the juncdb output dict of dicts
# dtrans - the dict of donors and acceptors. this is where
# we will count hits to those
def parse_sfish(f, djuncs, dtrans):

	szl = ""
	aln = []
	itpm = 0

	# open file
	fin = open(f, "r")

	# read header, find tpm column
	aln = fin.readline().strip().split("\t")
	while itpm < len(aln):
		if aln[itpm] == "TPM":
			break
		itpm += 1


	# read file
	for szl in fin:
		aln = szl.strip().split("\t")

		count = float(aln[itpm+1])

		if count > 0:

			if is_da(aln[0]):
				if aln[0] in dtrans:
					dtrans[aln[0]]['skip_count'] += count

			else:
				# not a donor/acceptor so this is a junction
				if aln[0] in djuncs:
					djuncs[aln[0]]['count'] += count

					if djuncs[aln[0]]['did'] in dtrans:
						dtrans[djuncs[aln[0]]['did']]['use_count'] += count

					if djuncs[aln[0]]['aid'] in dtrans:
						dtrans[djuncs[aln[0]]['aid']]['use_count'] += count


	fin.close()

	return 0


#
# Parse the juncdb gtf into a couple dicts. one is a dict by junction id
# where each element is a list containing the final output row for
# for output. the other is a table that translates the theta entries
# into their junction name buddies
def parse_gtf(fname, djuncs, dtrans):
	# variables
	gtfdb = {}
	grow = None

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:

		grow = GtfRow()
		grow.parse(szl.strip())

		if grow.db=="psi":

			# expand the tid out to get the actual splice site info
			r = re.search("^([^\:]+)\:([0-9]+)\-([0-9]+)$", grow.tid)
			jstart = int(r.group(2))
			jend = int(r.group(3))

			if jstart > jend:
				tmp = jstart
				jstart = jend
				jend = tmp

			# create donor and acceptor ids
			if grow.strand == "-":
				did = "{}:{}:3p".format(grow.rname, jstart)
				aid = "{}:{}:5p".format(grow.rname, jend)
			else:
				did = "{}:{}:5p".format(grow.rname, jstart)
				aid = "{}:{}:3p".format(grow.rname, jend)

			# this is a junction row so this is where we make an output list for 
			# the djuncs dict
			if grow.tid not in djuncs:
				djuncs[grow.tid] = dict(chrom=grow.rname, start=jstart, end=jend, 
					count=0.0, id="{}:{}-{}".format(grow.rname, jstart, jend), strand=grow.strand, donor_ovl=0.0, donor_use=0.0, 
					acc_ovl=0.0, acc_use=0.0, psi5="NA", psi3="NA", theta5="NA", theta3="NA", 
					gene_name=grow.attrs["gene_name"], gene_id=grow.attrs["gene_id"], 
					transcript_id=grow.attrs["oId"], did=did, aid=aid, jid=grow.tid)

			if did not in dtrans:
				dtrans[did] = dict(jid=grow.tid, use_count=0.0, skip_count=0.0)

			if aid not in dtrans:
				dtrans[aid] = dict(jid=grow.tid, use_count=0.0, skip_count=0.0)

	fin.close()

	return 0

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

#==============================================================================
#
# class definitions
#
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


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="This file is meant to run on the output of Sailfish when quantifying psi/theta junction database hits from raw reads.")
parser.add_argument('jdb_gtf', type=str, 
	help="Juncdb GTF annotation")
parser.add_argument('sfish_quant', type=str, 
	help="Sailfish quantification output.")

args = parser.parse_args()

#
# check args
if not isfile(args.jdb_gtf):
	sys.stderr.write("Error: input GTF file doesn't exist: {}\n".format(args.jdb_gtf))
	sys.exit(1)

if not isfile(args.sfish_quant):
	sys.stderr.write("Error: input Sailfish output doesn't exist: {}\n".format(args.sfish_quant))
	sys.exit(1)

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

