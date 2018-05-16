#!/usr/bin/python
#==============================================================================
# template.py
#
# Shawn Driscoll
# 20160712
# 
# Use mygene package to interface 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
import mygene
mg = mygene.MyGeneInfo()

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
	
	rres = None
	eid = 0
	dgenes = {}
	gid = ""
	gname = ""
	dkegg = {}

	sys.stderr.write("parsing info table\n")
	fin = open(args.infile, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		gname = aln[4]
		gid = aln[2]

		if gname not in dgenes:
			dgenes[gname] = { 'gene_id':[], 'id':[], 'alias':[], 'type':[], 'swiss.prot':[], 'trembl':[], 'unigene':[] }
	
		dgenes[gname]["gene_id"] += [gid]

	fin.close()

	sys.stderr.write("fetching mygene info\n")

	lgenes = dgenes.keys()
	n = len(lgenes)
	i = 0
	j = 0
	istep = 1000

	while i < n:
		j = i + istep
		if j > n:
			j = n

		# slice gene list keeping only unique gene names
		gset = list(set(lgenes[i:j]))

		# query
		rres = mg.querymany(gset, scopes="symbol", fields="entrezgene,symbol,alias,type_of_gene,uniprot,unigene", species=args.species)

		# loop through result
		if len(rres) > 0:
			for k in range(len(rres)):

				if "query" in rres[k]:

					gname = rres[k]['query']

					if gname not in dgenes:
						continue

					if "entrezgene" in rres[k]:
						dgenes[gname]['id'] += [rres[k]['entrezgene']]	
					elif "_id" in rres[k]:
						dgenes[gname]['id'] += [rres[k]['_id']]


					if "alias" in rres[k]:
						if type(rres[k]['alias']) is list:
							dgenes[gname]['alias'] += rres[k]['alias']
						else:
							dgenes[gname]['alias'].append(rres[k]['alias'])

					if "unigene" in rres[k]:
						if type(rres[k]['unigene']) is list:
							dgenes[gname]['unigene'] += rres[k]['unigene']
						else:
							dgenes[gname]['unigene'].append(rres[k]['unigene'])

					if "type_of_gene" in rres[k]:
						if type(rres[k]['type_of_gene']) is list:
							dgenes[gname]['type'] += rres[k]['type_of_gene']
						else:
							dgenes[gname]['type'].append(rres[k]['type_of_gene'])

					if 'uniprot' in rres[k]:
						if 'Swiss-Prot' in rres[k]['uniprot']:
							if type(rres[k]['uniprot']['Swiss-Prot']) is list:
								dgenes[gname]['swiss.prot'] += rres[k]['uniprot']['Swiss-Prot']
							else:
								dgenes[gname]['swiss.prot'].append(rres[k]['uniprot']['Swiss-Prot'])

						if 'TrEMBL' in rres[k]['uniprot']:
							if type(rres[k]['uniprot']['TrEMBL']) is list:
								dgenes[gname]['trembl'] += rres[k]['uniprot']['TrEMBL']
							else:
								dgenes[gname]['trembl'].append(rres[k]['uniprot']['TrEMBL'])

		i = j

	# print results by adding new columns to the input

	fin = open(args.infile, "r")
	fout = open("{}.ext".format(args.infile), "w")

	fout.write("\t".join(["chrom", "db", "gene_id", "transcript_id", "gene_name", "locus", "unigene", "swiss_prot", "TrEMBL", "eeid", "type", "alias"]))
	fout.write("\n")

	for szl in fin:
		aln = szl.strip().split("\t")

		gname = aln[4]
		dgene = dgenes[gname]

		gid = "u"
		galias = "u"
		gtype = "u"
		gswiss = "u"
		gtrembl = "u"
		guni = "u"

		if len(dgene['id']) > 0:
			gid = ",".join(map(str, uniq(dgene['id'])))

		if len(dgene['alias']) > 0:
			galias = ",".join(uniq(dgene['alias']))

		if len(dgene['type']) > 0:
			gtype = ",".join(uniq(dgene['type']))

		if len(dgene['swiss.prot']) > 0:
			gswiss = ",".join(uniq(dgene['swiss.prot']))

		if len(dgene['trembl']) > 0:
			gtrembl = ",".join(uniq(dgene['trembl']))

		if len(dgene['unigene']) > 0:
			guni = ",".join(uniq(dgene['unigene']))

		aln += [guni, gswiss, gtrembl, gid, gtype, galias]

		fout.write("\t".join(aln))
		fout.write("\n")

	fin.close()
	fout.close()

	return 0


def uniq(v):

	if type(v) is list:
		vhat = list(set(v))

	else:
		vhat = list(v)

	return vhat




#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Extends gene annotation info from mygene. Writes a new file (not to STDOUT)")
parser.add_argument('infile', type=str, help="Input file")
parser.add_argument("-s", "--species", type=str, default="mouse", 
	help="Species is one of human, mouse or rat [mouse]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


