#!/usr/bin/python
#==============================================================================
# template.py
#
# Shawn Driscoll
# 20160712
# 
# Use mygene package to interface with mygene.info database to fetch 
# KEGG pathway annotation for genes. genes are expected to be in the 
# <species_file>.info format (i.e. gene id is in column 3 and gene name 
# is in column 5)
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
	dinterpro = {}

	gtoi = "genes2interpro.tsv"
	itoi = "interproAnnot.tsv"

	sys.stderr.write("parsing info table\n")
	fin = open(args.infile, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		gname = aln[4]
		gid = aln[2]

		if gname not in dgenes:
			dgenes[gname] = { 'gene_id':set(), 'entrez':set(), 'interpro.id':set() }
	
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
		rres = mg.querymany(gset, scopes="symbol", fields="entrezgene,symbol,interpro", species=args.species)

		# loop through result
		if len(rres) > 0:
			for k in range(len(rres)):
				if "interpro" in rres[k]:

					gname = rres[k]['query']

					if type(rres[k]['interpro']) is dict:
						
						iid = rres[k]['interpro']['id']
						idesc = rres[k]['interpro']['desc']
						isdesc = rres[k]['interpro']['short_desc']
						
						if iid not in dinterpro:
							dinterpro[iid] = [isdesc, idesc]

						dgenes[gname]['interpro.id'].update([iid])

					else:
						# loop through list of terms

						for m in rres[k]['interpro']:

							iid = m['id']
							idesc = m['desc']
							isdesc = m['short_desc']
							
							if iid not in dinterpro:
								dinterpro[iid] = [isdesc, idesc]

							dgenes[gname]['interpro.id'].update([iid])


		i = j



	# print results
	fout = open(gtoi, "w")
	# write header
	fout.write("gene_name\tgene_id\tid\n")

	for gname in dgenes.keys():
		if len(dgenes[gname]['interpro.id']) > 0:

			# print one line per gene - interpro id association like a graph table
			lkid = list(dgenes[gname]['interpro.id'])
			for k in lkid:
				lout = [gname, ",".join(list(dgenes[gname]['gene_id'])), k] 
				fout.write("\t".join(lout))
				fout.write("\n")


	fout.close()

	# print annotation to a separate file

	fout = open(itoi, "w")
	# write header
	fout.write("id\tshort_desc\tdesc\n")
	for iid in dinterpro.keys():
		lout = [iid] + dinterpro[iid]
		fout.write("\t".join(lout))
		fout.write("\n")

	fout.close()

	return 0





#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('infile', type=str, help="Input file")
parser.add_argument("-s", "--species", type=str, default="mouse", 
	help="Species name {human, mouse, rat}")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


