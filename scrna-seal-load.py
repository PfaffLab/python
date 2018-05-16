#!/usr/bin/python
#==============================================================================
# scrna-seal-load.py
#
# Shawn Driscoll
# 20170912
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse the output of several scrna-seal-count.py ouputs and load counts 
# into a table
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from os.path import isfile, expanduser, basename, dirname
from collections import defaultdict
from time import localtime, time
#from multiprocessing import Pool
from os import system
from hashlib import md5
import numpy as np

# from igraph import *
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

np.random.seed(123)

HOME = expanduser("~")


#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	bc_list = []
	genes = []
	gene_index = defaultdict(int)
	counts = []
	err = []
	num_files = len(args.umi_counts)
	has_err = False
	any_has_err = False
	cell_sums = [0 for i in range(num_files)]

	super_total_counts = 0

	cell_file = "scrna_seal_cells.tsv"
	gene_file = "scrna_seal_genes.tsv"
	counts_file = "scrna_seal_countsMM.tsv"
	err_file = "scrna_seal_errEstMM.tsv"

	# 
	# first let's pick apart the file names since they include the cell barcode
	bc_list = ["" for i in range(num_files)]
	i = 0
	for f in args.umi_counts:
		r = re.search("\_([ACTG]+)\.umi_counts$", f)
		if r:
			# name it a-la cellranger
			bc_list[i] = r.group(1) + "-1"
		else:
			# if we can't extract the barcode then just use the whole file name
			# without extension as the cell id
			bc_list[i] = basename(f).split(".")[0]

		i += 1

	t0 = time()
	message("Writing cells file")
	szout = ""
	for k in bc_list:
		szout += "{}\n".format(k)
	with open(cell_file, "w") as fout:
		fout.write(szout)
	sys.stderr.write("{} sec\n".format(time()-t0))

	progress_message("Loading {}".format(args.umi_counts[0]))

	lcounts = []
	lerr = []

	with open(args.umi_counts[0], "r") as fin:

		# header. figure out if the file contains resampling uncertainty estimates
		szl = fin.readline()
		aln = szl.strip().split("\t")
		if len(aln) > 3:
			# this file has resampling uncertainty
			has_err = True
			any_has_err = True

		i = 0
		for szl in fin:
			aln = szl.strip().split("\t")
			gene_index[aln[0]] = i
			genes.append(aln[0])
			idx = i
			count = float(aln[2])

			i += 1

			if count > 0:
				super_total_counts += count
				# keep this line
				lcounts.append([idx, 0, count])
				if has_err:
					lerr.append([idx, 0, aln[4]])

	# total count of genes discovered
	num_genes = len(genes)

	# loop through remaining files
	for i in range(1, num_files):
		has_err = False
		
		progress_message("Loading {}".format(args.umi_counts[i]))
		with open(args.umi_counts[i], "r") as fin:
			szl = fin.readline()
			aln = szl.strip().split("\t")
			if len(aln) > 3:
				# this file has resampling uncertainty
				err[i] = [0 for j in range(num_genes)]
				has_err = True
				any_has_err = True

			for szl in fin:
				aln = szl.strip().split("\t")
				idx = gene_index[aln[0]]
				count = float(aln[2])
				if count > 0:
					super_total_counts += count
					lcounts.append([idx, i, count])
					if has_err:
						lerr.append([idx, i, float(aln[4])])

	# now we have the three column counts table. first index is the cell, second is the gene
	# we can calculate the normalization factors.

	sys.stderr.write("\n")
	sys.stderr.write("{} sec\n".format(time()-t0))

	if False:

		temp = sorted(cell_sums)
		if num_files % 2 == 0:
			med_count = (cell_sums[num_files/2] + cell_sums[num_files/2+1])*1.0/2
		else:
			med_count = cell_sums[int(math.floor(num_files*1.0/2)+1)]
		
		sizes = [a/med_count for a in cell_sums]

		norm_counts = [[] for i in range(num_files)]
		for i in range(num_files):
			norm_counts[i] = [x/sizes[i] for x in counts[i]]

	t0 = time()
	message("Building output table")

	# make output
	szout = "{}\t{}\t{}\n".format(num_genes, num_files, super_total_counts)
	for l in lcounts:
		szout += "\t".join(map(str, l))
		szout += "\n"

	with open(counts_file, "w") as fout:
		fout.write(szout)

	# write genes file
	szout = ""
	for g in genes:
		szout += "{}\n".format(g)

	with open(gene_file, "w") as fout:
		fout.write(szout)

	# do we have error estimates?
	if len(lerr) > 0:
		szout = "{}\t{}\t{}\n".format(num_genes, num_files, 0)
		for l in lerr:
			szout += "\t".join(map(str, l))
			szout += "\n"

		with open(err_file, "w") as fout:
			fout.write(szout)

	sys.stderr.write("{} sec\n".format(time()-t0))

	return 0
		

def quick_gtf_parse(f):
	
	attrs = {}
	annot = {}
	gname2tid = {}

	with open(f, "r") as fin:

		for szl in fin:
			
			# extract transcript id 
			r = re.search("transcript_id \"([^\"]+)\"", szl)
			if r:
				tid = r.group(1)
			else:
				sys.stderr.write("Warning: failed to parse transcript id @ {}".format(szl))

			# if transcript id is already in the annot dict then we don't have to do anything 
			# further with this line
			if tid not in annot:
				# split the line up and parse out attributes to extract gene name and gene id
				aln = szl.strip().split("\t")

				# initalize attributes dict
				attrs = {}

				temp = aln[8].split(";")
				for k in temp:
					# remove whitespace
					k = re.sub("^[\s]+", "", k.strip())
					# continue if we have characters left
					
					pos = k.find(" ")
					if len(k) > 3 and pos > 0:
						# copy out the attribute name
						aid = k[0:pos]
						# copy out the attribute value
						aval = k[(pos+1):len(k)]
						# remove quotes if they are there
						if aval[0]=="\"":
							aval = aval[1:(len(aval)-1)]
						
						attrs[aid] = aval
				
				# add to the annotation dict
				annot[tid] = {
					'chrom': aln[0], 
					'db': aln[1],
					'feature': aln[2], 
					'strand': aln[6],
					'transcript_id': tid, 
					'gene_id': attrs['gene_id']
				}
				
				if 'gene_name' in attrs:
					annot[tid]['gene_name'] = attrs['gene_name']
				else:
					annot[tid]['gene_name'] = attrs['gene_id']
				
				# add translation from gene name back to transcript id
				gname2tid[annot[tid]["gene_name"]] = tid
				

	return annot, gname2tid

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0

def message(sz, show_date=True):
	if show_date:
		msg = "[{}] {}".format(time_string(), sz)
	else:
		msg = sz

	sys.stderr.write(sz + "\n")

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def print_exception():
	exc_type, exc_value, exc_traceback = sys.exc_info()
	print "*** print_tb:"
	traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
	print "*** print_exception:"
	traceback.print_exception(exc_type, exc_value, exc_traceback,
	                          limit=2, file=sys.stdout)
	print "*** print_exc:"
	traceback.print_exc()
	print "*** format_exc, first and last line:"
	formatted_lines = traceback.format_exc().splitlines()
	print formatted_lines[0]
	print formatted_lines[-1]
	print "*** format_exception:"
	print repr(traceback.format_exception(exc_type, exc_value,
	                                      exc_traceback))
	print "*** extract_tb:"
	print repr(traceback.extract_tb(exc_traceback))
	print "*** format_tb:"
	print repr(traceback.format_tb(exc_traceback))
	print "*** tb_lineno:", exc_traceback.tb_lineno	

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Load the results from scrna-seal-count.py runs into a single table.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('umi_counts', nargs="+", 
	help="Two or more outputs from scrna-seal-count.py")

annot = parser.add_mutually_exclusive_group(required=False)
annot.add_argument('-g', type=str, help="GTF annotation to expand on gene annotation")
annot.add_argument('-i', type=str, help="Parsed GTF 'info' file to use to expand annotation of genes")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

