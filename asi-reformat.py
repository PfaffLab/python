#!/usr/bin/env python
#==============================================================================
# asi-reformat.py
#
# Shawn Driscoll
# 20160715
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Reformat the ASI quantifications from asi-count-from-junctions.py which 
# bundles locus type groups of alternative paths rather than paired off 
# sets. 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
import copy

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
	daloc = {}
	dheader = {}
	lheader = []
	aeidBase = "AEID_"
	iretBase = "IRET_"
	aeidx = 0

	if not isfile(args.asi):
		sys.stderr.write("Input ASI file does not exist ({})".format(args.asi))
		return 1

	sys.stderr.write("Parsing input file ({})\n".format(args.asi))

	# open file and parse it into the dict
	fin = open(args.asi, "r")

	# read in header line
	lheader = fin.readline().strip().split("\t")
	for i in range(len(lheader)):
		dheader[lheader[i]] = i

	for szl in fin:
		aln = szl.strip().split("\t")
		daln = parse_line_to_dict(aln, dheader)

		aloc = daln['aloc_id']
		if aloc not in daloc:
			daloc[aloc] = {'type': daln["asi_type"], 
				'gene_name': daln["gene_name"], 
				'gene_id': daln["gene_id"], 'paths': []}

		daloc[aloc]['paths'].append(daln)

	fin.close()

	# --
	# parsing of the file is complete. now each aloc has to be transformed into 
	# single line even rows.
	# --

	sys.stderr.write("Pairing alternative splice events and printing to STDOUT\n")

	for aloc in sorted(daloc.keys()):

		# number of features
		n = len(daloc[aloc]['paths'])
		
		# make pairwise merges
		if daloc[aloc]['type'] == "alt.cassette":
			for i in range(1, n):
				aeidx += 1
				aeid = "{}{:08d}".format(aeidBase, aeidx)

				aeout = merge_cassette(daloc[aloc]['paths'][0], daloc[aloc]['paths'][i])
				aeout.aeid = aeid

				if aeidx==1:
					lout = aeout.header_row()
					print "\t".join(lout)

				lout = aeout.tolist()
				print "\t".join(map(str, lout))

		else:
			for i in range(0, n-1):
				for j in range(i+1, n):
					aeidx += 1

					aeid = "{}{:08d}".format(aeidBase, aeidx)
					aetype = daloc[aloc]['type']
					
					if aetype == "alt.3p" or aetype == "alt.5p":
						aeout = merge_altp(daloc[aloc]['paths'][i], daloc[aloc]['paths'][j])
					elif aetype == "alt.end" or aetype == "alt.start":
						aeout = merge_alt_start_end(daloc[aloc]['paths'][i], daloc[aloc]['paths'][j])
					elif aetype == "mut.ex":
						aeout = merge_mutex(daloc[aloc]['paths'][i], daloc[aloc]['paths'][j])
					
					aeout.aeid = aeid

					if aeidx==1:
						lout = aeout.header_row()
						print "\t".join(lout)

					lout = aeout.tolist()
					print "\t".join(map(str, lout))


	if args.juncdb is not None:
		if isfile(args.juncdb):

			ljparse = []

			# parse the juncdb file and make intron.retention lines
			#sys.stderr.write("Parsing juncdb file ({}) to append intron retention events\n".format(args.juncdb))

			fin = open(args.juncdb, "r")
			# pull in header
			lheader = fin.readline().strip().split("\t")
			dheader = {}
			for i in range(len(lheader)):
				dheader[lheader[i]] = i

			# parse in the file and sort the junctions by their id so we get a kind of deterministic 
			# output
			sys.stderr.write("Parsing {}\n".format(args.juncdb))
			for szl in fin:
				ljparse.append(parse_line_to_dict(szl.strip().split("\t"), dheader))

			fin.close()

			# sorting...

			j = 0
#			for szl in fin:
			for daln in sorted(ljparse, key=lambda x: x['id']):
				#daln = parse_line_to_dict(szl.strip().split("\t"), dheader)

				aeidx += 1
				j += 1 # count introns

				aeid = "{}{:08d}".format(aeidBase, aeidx)
				aeout = AERow(aeid=aeid, 
					alocid="{}{:08d}".format(iretBase, j),
					gid=daln['gene_id'], 
					gname=daln['gene_name'], 
					aetype='intron.retention',
					strand=daln['strand'],
					intronsA=daln['id'],
					intronsB="None",
					tidA=daln['transcript_id'],
					tidB="None",
					countA=daln['count'],
					countB=float(float(daln['donor_ovl'])+float(daln['acc_ovl']))/2.0)

				lout = aeout.tolist()

				print "\t".join(map(str, lout))

	sys.stderr.write("Done.\n")

	return 0

#
# a and b are dicts and would be either an alt 3' or 5' feature. for condition 
# a we want to use the one with the shorter intron
def merge_altp(a, b):

	ilen1 = length_from_position(a['introns'])
	ilen2 = length_from_position(b['introns'])

	if ilen2 < ilen1:
		ref = b
		alt = a
	else:
		ref = a
		alt = b

	aeout = AERow(alocid=ref['aloc_id'], 
		gid=ref['gene_id'], 
		gname=ref['gene_name'], 
		aetype=ref['asi_type'], 
		strand=ref['strand'], 
		intronsA=ref['introns'], 
		intronsB=alt['introns'], 
		tidA=ref['transcript_ids'], 
		tidB=alt['transcript_ids'], 
		countA=ref['mean_path_count'], 
		countB=alt['mean_path_count'])

	return aeout

# 
# feature with the longer intron is reference
def merge_alt_start_end(a, b):

	ilen1 = length_from_position(a['introns'])
	ilen2 = length_from_position(b['introns'])

	if ilen2 < ilen1:
		ref = a
		alt = b
	else:
		ref = b
		alt = a

	aeout = AERow(alocid=ref['aloc_id'], 
		gid=ref['gene_id'], 
		gname=ref['gene_name'], 
		aetype=ref['asi_type'], 
		strand=ref['strand'], 
		intronsA=ref['introns'], 
		intronsB=alt['introns'], 
		tidA=ref['transcript_ids'], 
		tidB=alt['transcript_ids'], 
		countA=ref['mean_path_count'], 
		countB=alt['mean_path_count'])

	return aeout


# 
# feature with no exon (i.e. one intron) is reference.
def merge_cassette(a, b):

	num_introns_a = len(a['introns'].split("|"))
	num_introns_b = len(b['introns'].split("|"))

	if num_introns_a==1:
		ref = a
		alt = b
		alt_count = num_introns_b
	else:
		ref = b
		alt = a
		alt_count = num_introns_a

	aetype = "alt.cassette"
	if alt_count > 2:
		aetype = "tandem.cassette"

	aeout = AERow(alocid=ref['aloc_id'], 
		gid=ref['gene_id'], 
		gname=ref['gene_name'], 
		aetype=aetype, 
		strand=ref['strand'], 
		intronsA=ref['introns'], 
		intronsB=alt['introns'], 
		tidA=ref['transcript_ids'], 
		tidB=alt['transcript_ids'], 
		countA=ref['mean_path_count'], 
		countB=alt['mean_path_count'])

	return aeout


# 
# merge them as they come
def merge_mutex(a, b):

	ref = a
	alt = b

	aeout = AERow(alocid=ref['aloc_id'], 
		gid=ref['gene_id'], 
		gname=ref['gene_name'], 
		aetype=ref['asi_type'], 
		strand=ref['strand'], 
		intronsA=ref['introns'], 
		intronsB=alt['introns'], 
		tidA=ref['transcript_ids'], 
		tidB=alt['transcript_ids'], 
		countA=ref['mean_path_count'], 
		countB=alt['mean_path_count'])

	return aeout

def parse_line_to_dict(ll, dheader):
	dout = {}

	for k in dheader.keys():
		dout[k] = ll[dheader[k]]

	return dout

# --
# calculate length from the position string
def length_from_position(sz):
	ll = explode_position(sz)
	return ll[2]-ll[1]

# ---
# explode position string to chrom/start/end
def explode_position(sz):
	lout = ["", 0, 0]
	ls1 = sz.split(":")
	ls2 = ls1[1].split("-")
	lout[0] = ls1[0]
	lout[1] = int(ls2[0])
	lout[2] = int(ls2[1])
	return lout


# --
# class for the merged output rows
# --

class AERow(object):
	def __init__(self, aeid="", alocid="", gid="", gname="", aetype="", strand="", intronsA="", 
		intronsB="", tidA="", tidB="", countA=0, countB=0, psA=0.0):

		self.aeid = aeid
		self.alocid = alocid
		self.gid = gid
		self.gname = gname
		self.aetype = aetype
		self.strand = strand
		self.intronsA = intronsA
		self.intronsB = intronsB
		self.tidA = tidA
		self.tidB = tidB
		self.countA = float(countA)
		self.countB = float(countB)
		self.psA = psA

		return None

	def header_row(self):

		lout = ["aeid", 
			"aloc_id", "gene_id", "gene_name", "ae_type", "strand", 
			"introns_a", "transcript_id_a", "introns_b", "transcript_id_b", 
			"count_a", "count_b", "total_count", "ps_a"]

		return lout		

	# 
	# this function parses the tracking file row into this class
	def tolist(self):

		totalCount = self.countA+self.countB
		if totalCount > 0:
			self.psA = self.countA/totalCount
		else:
			self.psA = -1

		lout = [
			self.aeid, 
			self.alocid, 
			self.gid, 
			self.gname, 
			self.aetype, 
			self.strand,
			self.intronsA,
			self.tidA,
			self.intronsB,
			self.tidB,
			self.countA,
			self.countB,
			totalCount,
			self.psA
		]

		return lout


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('asi', type=str, help="asi quantification file from asi-count-from-junctions.py")
parser.add_argument('-j', '--juncdb', type=str, default=None, 
	help="Juncdb quantification. This will merge intron retention into the output.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

