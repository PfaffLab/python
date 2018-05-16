#!/usr/bin/python
#==============================================================================
# annotate-kallisto-pseudo.py
#
# Shawn Driscoll
# 20171205
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Brainstorm:
# kallisto output is maybe easiest to load in R but R is painfully slow at it
# so here's what we have to do:
# load the FAI file, parse the transcript names from the first column into a 
# list. The "matrix.ec" file will reference these names by zero-based index.
# load an annotation file that translates transcript ids into gene names. 
#==============================================================================

import sys
import argparse
import re
import traceback
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime, time
from math import floor
from random import random
#import subprocess as sp
#from os import system
#from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock

# from igraph import *
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

ANNOT_GID = 2
ANNOT_TID = 3
ANNOT_GENE_NAME = 4
ANNOT_GENE_LOCUS = 5

#==============================================================================
# main
#==============================================================================

def main(args):


	# t0 = time()
	#
	# let user know what's up
	# sys.stderr.write("{} sec\n".format(time()-t0))

	tcc = "{}/matrix.tsv".format(args.source_path)
	ec = "{}/matrix.ec".format(args.source_path)
	cells = "{}/matrix.cells".format(args.source_path)
	gene_to_ec = defaultdict(list)
	gene_master_list = []
	gene_master_lk = {}
	gene_master_idx = 0

	##
	## check all files
	##
	file_missing_flag = False

	if not isfile(args.fai):
		file_missing_error(args.fai)
		file_missing_flag = True

	if not isfile(args.annotation):
		file_missing_error(args.fai)
		file_missing_flag = True

	if not isfile("{}/matrix.tsv".format(args.source_path)):
		file_missing_error("{}/matrix.tsv".format(args.source_path))
		file_missing_flag = True

	if not isfile("{}/matrix.ec".format(args.source_path)):
		file_missing_error("{}/matrix.ec".format(args.source_path))
		file_missing_flag = True

	if not isfile("{}/matrix.cells".format(args.source_path)):
		file_missing_error("{}/matrix.cells".format(args.source_path))
		file_missing_flag = True

	if file_missing_flag:
		return 1

	##
	## load FAI
	sys.stderr.write("Loading FAI...\n")
	t0 = time()
	fai = parse_fai(args.fai)
	sys.stderr.write("{} sec\n".format(time()-t0))

	##
	## load annotation
	sys.stderr.write("Loading annotation...\n")
	t0 = time()
	dannot = parse_annotation(args.annotation)
	sys.stderr.write("{} sec\n".format(time()-t0))

	##
	## load cells file
	##
	sys.stderr.write("Loading cells file...\n")
	t0 = time()
	cells_list = []
	with open(cells, "r") as fin:
		for szl in fin:
			cells_list.append(szl.strip())

	num_cells = len(cells_list)

	sys.stderr.write("{} sec\n".format(time()-t0))
	sys.stderr.write("number of cells:     {}\n".format(num_cells))

	##
	## load ec file
	sys.stderr.write("Loading equivalence class matrix and annotating...\n")
	t0 = time()
	
	ec_tid = []
	ec_gname = []

	# list of ec indices that represent only a single gene
	ec_single_gene = []
	ec_multi_gene = []

	with open(ec, "r") as fin:
		for szl in fin:

			aln = szl.strip().split("\t")
			ec_idx = int(aln[0])
			idx = map(int, aln[1].split(","))
			ec_tid.append([fai[i] for i in idx])
			
			gnames = []
			for tid in ec_tid[-1]:
				if tid in dannot:
					gnames.append(dannot[tid][ANNOT_GENE_NAME])
			
			tmp = list(set(gnames))
			if len(tmp)==1:
				ec_single_gene.append(ec_idx)
				gene_to_ec[tmp[0]].append(ec_idx)

				if tmp[0] not in gene_master_lk:
					gene_master_lk[tmp[0]] = gene_master_idx
					gene_master_list.append(tmp[0])
					gene_master_idx += 1

			else:
				ec_multi_gene.append(ec_idx)

			ec_gname.append(list(set(gnames)))

	sys.stderr.write("{} sec\n".format(time()-t0))

	num_ec = len(ec_tid)
	num_genes = len(gene_master_list)

	sys.stderr.write("equivalence classes: {}\n".format(len(ec_tid)))
	sys.stderr.write("single gene ec's:    {}\n".format(len(ec_single_gene)))

	#
	##
	### now we can parse the tcc matrix and figure this junk out
	##
	#

	sys.stderr.write("Loading ec umi count matrix...\n")
	t0 = time()

	# setup lists to gather umi counts:
	# ec level
	umi = [[] for i in range(num_cells)]
	# gene level	
	gene_umi = [[] for i in range(num_cells)]
	cell_umi_stats = [[0, 0] for i in range(num_cells)]

	with open(tcc, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			
			ec_id = int(aln[0])
			cell_id = int(aln[1])
			umi_count = int(aln[2])

			# columns: ec_id  cell_id  umi_count
			if len(umi[cell_id]) == 0:
				umi[cell_id] = [0 for i in range(num_ec)]
				gene_umi[cell_id] = [0 for i in range(num_genes)]

			umi[cell_id][ec_id] += umi_count
			# if this ec is a single gene ec then update the gene's count
			if len(ec_gname[ec_id])==1:
				gname = ec_gname[ec_id][0]
				gene_umi[cell_id][gene_master_lk[gname]] += umi_count

			# update total umi count for this cell
			cell_umi_stats[cell_id][0] += umi_count

	sys.stderr.write("{} sec\n".format(time()-t0))

	##
	## write the gene level umi counts to a file. these are umi counts discarding any
	## potentially ambiguous hits
	##

	sys.stderr.write("Writing gene level UMI counts to file...\n")
	t0 = time()

	szout = "id\t" + "\t".join(cells_list) + "\n"
	for i in range(num_genes):
		szout += gene_master_list[i]
		for cellid in range(num_cells):
			szout += "\t{:d}".format(gene_umi[cellid][i])
		szout += "\n"

	with open("{}/gene_level_noambig.tsv".format(args.source_path), "w") as fout:
		fout.write(szout)

	sys.stderr.write("{} sec\n".format(time()-t0))



	##
	## now we are gonna try something funny.  let's check out those EC's with multiple
	## genes that have hits. for each one pull the hits to those genes from ec's with 
	## single genes.  establish a weighting based on those hits in order to divide up 
	## the shared hits.
	##

	sys.stderr.write("Initializing weights for ambiguous count assignment...\n")
	t0 = time()

	# copy original and keep it
	gene_umi_base = copy_expression_matrix(gene_umi)
	
	# in this copy we will update the counts at each iteration. the base
	# 'gene_umi' will be used to calculate weights. at the end of each loop
	# we have to copy gene_umi0 to gene_umi and then copy gene_umi_base
	# to gene_umi0
	gene_umi0 = copy_expression_matrix(gene_umi)
	

	# create initial weights for multi-target EC. initial weight for each gene is 1/n
	# where n = number of genes in the EC
	mweights = [[] for i in range(num_cells)]
	for cellid in range(num_cells):
		progress_message("initalizing cell {}".format(cellid))
		mweights[cellid] = [None for i in range(num_ec)]
		for ecid in ec_multi_gene:
			if umi[cellid][ecid] > 0:
				n = len(ec_gname[ecid])
				mweights[cellid][ecid] = [1.0/n for i in range(n)]
	
	sys.stderr.write("\n")
	weightse = [0 for i in range(num_ec)]
	
	sys.stderr.write("{} sec\n".format(time()-t0))
	sys.stderr.write("Assigning ambiguous UMI counts...\n")
	t0 = time()
	
	max_iter = 100
	tol = 1e-3
	iter = 0
	
	#watchlist = [380869, 384201]
	
	while True:
		sosq = 0
		iter += 1
		
		##
		## note that as this loop proceeds the weight vector at a particular 
		## ecid would change as we modify the base counts. so we need to reference
		## the copy of the base counts for each update
		for cellid in range(num_cells):
			progress_message("processing cell {}".format(cellid))
			for ecid in ec_multi_gene:
				# check umi count for this multi-gene hit
				if umi[cellid][ecid] > 0:
					# positive umi, get the gene level counts which we'll use as weights to 
					# partition the shared hits.
					gcounts = []
					for gname in ec_gname[ecid]:
						gcounts.append(gene_umi[cellid][gene_master_lk[gname]])
	
					scount = sum(gcounts)
					if scount > 0:
						# we have counts. make the weights and then add count into 
						# the gene-level counts at each gene
						w = [a*1.0/scount for a in gcounts]
						# calculate change from previous pass
						w0 = mweights[cellid][ecid]
						#if random() < 0.01:
						#if ecid in watchlist:
						#	print cellid, ecid, w0, w

						sosq += sum([(w[i]-w0[i])**2 for i in range(len(w))])
						mweights[cellid][ecid] = list(w)
						
						# update counts
						for k in range(len(w)):
							gname = ec_gname[ecid][k]
							gidx = gene_master_lk[gname]
							gene_umi0[cellid][gidx] += w[k]*umi[cellid][ecid]
		
		# here we copy the updated counts to gene_umi and copy the base counts to 
		# gene_umi0
		gene_umi = copy_expression_matrix(gene_umi0)
		gene_umi0 = copy_expression_matrix(gene_umi_base)
		
		sosq = sosq*1.0/num_cells
		sys.stderr.write("\n")
		sys.stderr.write("Iterations: {}; sse: {}\n".format(iter, sosq))
		
		if iter >= max_iter or sosq < tol:
			break

	sys.stderr.write("{} sec\n".format(time()-t0))

	sys.stderr.write("Writing adjusted gene level UMI counts to file...\n")
	t0 = time()

	szout = "id\t" + "\t".join(cells_list) + "\n"
	for i in range(num_genes):
		szout += gene_master_list[i]
		for cellid in range(num_cells):
			szout += "\t{:d}".format(int(round(gene_umi[cellid][i])))
		szout += "\n"

	with open("{}/gene_level.tsv".format(args.source_path), "w") as fout:
		fout.write(szout)

	sys.stderr.write("{} sec\n".format(time()-t0))

	return 0

def copy_expression_matrix(m):

	ncol = len(m)
	nrow = len(m[0])

	lout = [[] for i in range(ncol)]
	for i in range(ncol):
		lout[i] = [m[i][j] for j in range(nrow)]

	return lout


def round(x):

	xhat = floor(x)

	if x - xhat >= 0.5:
		xhat += 1

	return xhat


def parse_fai(fname):

	tid = []
	with open(fname, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			tid.append(aln[0])

	return tid

def parse_annotation(fname):

	dannot = defaultdict(list)

	with open(fname, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			if aln[0] == "chrom":
				continue

			dannot[aln[ANNOT_TID]] = list(aln)

	return dannot

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0

def file_missing_error(fname):
	sys.stderr.write("Error: file {} is missing\n".format(fname))

def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

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

parser = argparse.ArgumentParser(description="Annotate the kallisto-pseudo 'batch' output for downstream analysis.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('source_path', type=str, help="Path to kallisto-pseudo output files")
parser.add_argument("fai", type=str, help="FAI index for FASTA that the kallisto index was built from")
parser.add_argument("annotation", type=str, help="Tab-delim file generated from 'parse-gtf' typically used as annotation")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

