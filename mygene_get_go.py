#!/usr/bin/python
#==============================================================================
# mygene_get_go.py
#
# Shawn Driscoll
# 20160712
# 
# Use mygene package to interface with mygene.info database to fetch 
# gene ontology (GO) for all genes in one of my ".info" files which are
# generated from GTF files
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
	dannot = {}

	gtoid_out = "genes2go.tsv"
	annot_out = "goAnnot.tsv"

	sys.stderr.write("parsing info table\n")
	fin = open(args.infile, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		gname = aln[4]
		gid = aln[2]

		if gname not in dgenes:
			dgenes[gname] = { 'gene_id':set(), 'id':set() }
	
		dgenes[gname]["gene_id"].update([gid])

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

		# slice gene list
		gset = lgenes[i:j]

		# query
		rres = mg.querymany(gset, scopes="symbol", fields="go", species=args.species)

		# loop through result
		if len(rres) > 0:
			for k in range(len(rres)):
				if "go" in rres[k]:

					gname = rres[k]['query']

					for go_cat in rres[k]['go'].keys():


						if type(rres[k]['go'][go_cat]) is dict:
							tid = rres[k]['go'][go_cat]['id']
							tname = rres[k]['go'][go_cat]['term']

							if tid not in dannot:
								dannot[tid] = [go_cat, tname]

							dgenes[gname]['id'].update([tid])

						else:
							# loop through list of pathways

							for m in rres[k]['go'][go_cat]:
								tid = m['id']
								tname = m['term']
								if tid not in dannot:
									dannot[tid] = [go_cat, tname]

								dgenes[gname]['id'].update([tid])


		i = j

	if False:

		for gname in dgenes.keys():
			# first query to get the entrez id or ensembl id for this gene 
			# symbol
			rres = mg.query("symbol:{}".format(gname), species=args.species)

			if 'hits' in rres:
				if len(rres['hits']) > 0:
					# we have a hit!
					eid = rres['hits'][0]['_id']
					# get desired annotation
					rres2 = mg.getgene(eid, "symbol,pathway.kegg")

					if rres2 is None:
						continue

					# if there was kegg info in there...
					if "pathway.kegg" in rres2:

						# if there was only a single pathway then it comes back
						# as a dict instead of a list
						if type(rres2['pathway.kegg']) is dict:
							kid = rres2['pathway.kegg']['id']
							kname = rres2['pathway.kegg']['name']

							if kid not in dkegg:
								dkegg[kid] = kname

							dgenes[gname]['kegg.id'].update([kid])

						else:
							# loop through list of pathways

							for k in rres2['pathway.kegg']:
								kid = k['id']
								kname = k['name']
								if kid not in dkegg:
									dkegg[kid] = kname

								dgenes[gname]['kegg.id'].update([kid])

	# print results
	fout = open(gtoid_out, "w")
	fout.write("gene_name\tgene_id\tid\n")
	for gname in dgenes.keys():
		if len(dgenes[gname]['id']) > 0:

			# print one line per gene - kegg id association like a graph table
			lkid = list(dgenes[gname]['id'])
			for k in lkid:
				lout = [gname, ",".join(list(dgenes[gname]['gene_id'])), k] 
				fout.write("\t".join(lout))
				fout.write("\n")


	fout.close()

	# print annotation to a separate file

	fout = open(annot_out, "w")
	fout.write("id\tnamespace\tname\n")
	for tid in dannot.keys():
		lout = [tid, "\t".join(dannot[tid])]
		fout.write("\t".join(lout))
		fout.write("\n")

	fout.close()

	return 0





#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('infile', type=str, help="Input file")
parser.add_argument('-s', '--species', default="mouse", type=str, 
	help="Species to fetch.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


