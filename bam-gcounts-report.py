#!/usr/bin/python
#==============================================================================
# scrna-bd-combine-results.py
#
# Shawn Driscoll
# 20170615
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Combine quantification results from multiple scrna-gcounts-bd runs
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime

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

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	n = len(args.quants)
	allQ = QuantCatalog()
	
	if n==1:
		error_message("This is only useful if you have more than one file!")
		return 1
	
	stub = "quant_"

	message("Loading files")
	for f in args.quants:
		progress_message("{}".format(f))
		quant = Quant()
		quant.parse(f)
		quant.summarize_gene_id_level()
		allQ.add_quant(quant)
	
	sys.stderr.write("\n")
	
	
	# all quants are loaded. now make some stats and generate the quantificaiton 
	# summary tables
	message("building tables")
	
	allQ.build_tables("gene_id")
	
	#
	# write tables
	#
	
	message("Writing Raw Counts Table")
	with open("{}raw_counts.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.counts)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.counts[i])))
			fout.write("\n")

	message("Writing Effective Counts Table")
	with open("{}effective_counts.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.eff_counts)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.eff_counts[i])))
			fout.write("\n")

	message("Writing FPKM Table")
	with open("{}FPKM.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.fpkm)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.fpkm[i])))
			fout.write("\n")

	message("Writing TPM Table")
	with open("{}TPM.tsv".format(stub), "w") as fout:
		
		fout.write("\t".join(allQ.annot_header + allQ.data_header) + "\n")
		for i in range(len(allQ.tpm)):
			fout.write("\t".join(allQ.annot[i] + map(str, allQ.tpm[i])))
			fout.write("\n")
	
#	message("Writing summary report")
#	
#	rpt = allQ.build_report()
#	with open("{}report.tsv".format(stub), "w") as fout:
#		for l in rpt:
#			fout.write("\t".join(map(str, l)) + "\n")

	return 0

class QuantCatalog(object):
	def __init__(self):
		self.quants = []
		self.num_quants = 0
	
		self.annot = None
		self.annot_header = None
		
		self.counts = None
		self.eff_counts = None
		self.fpkm = None
		self.tpm = None
		
		self.data_header = None
		
		
	def add_quant(self, obj):
		self.num_quants += 1
		self.quants.append(obj)
	
	def build_tables(self, type):
		if type=="gene":
			self.build_gene_level_tables()
		elif type=="gene_id":
			self.build_gene_id_level_tables()
	
	def build_gene_level_tables(self):
		
		if self.num_quants < 1:
			error_message("No samples!")
			return None		
	
		# build union of gene names from each sample
		gnames = set()
		for obj in self.quants:
			if obj.gn2idx is not None:
				gnames.update(obj.gn2idx.keys())
		
		gnames = sorted(list(gnames))
		num_genes = len(gnames)
		
		# now build the output
		tid = [None for i in range(num_genes)]		
		gid = [None for i in range(num_genes)]
		
		self.counts = [[0 for i in range(self.num_quants)] for j in range(num_genes)]
		self.eff_counts = [[0 for i in range(self.num_quants)] for j in range(num_genes)]
		self.fpkm = [[0 for i in range(self.num_quants)] for j in range(num_genes)]
		self.tpm = [[0 for i in range(self.num_quants)] for j in range(num_genes)]
				
		for i in range(num_genes):
			gname = gnames[i]
			for j in range(self.num_quants):
				obj = self.quants[j]
				if gname in obj.gn2idx:
					gidx = obj.gn2idx[gname]
						
					tid[i] = ";".join(obj.gene_annot[gname]['tid'])
					gid[i] = ";".join(obj.gene_annot[gname]['gid'])
					
					self.counts[i][j] = obj.gene_hits[gidx]
					self.eff_counts[i][j] = obj.gene_eff_hits[gidx]
					self.fpkm[i][j] = obj.gene_fpkm[gidx]
					self.tpm[i][j] = obj.gene_tpm[gidx]
		
		self.annot = [[] for i in range(num_genes)]
		self.annot_header = ["gene_name", "gene_id", "transcript_id"]
		for i in range(num_genes):
			self.annot[i] = [gnames[i], gid[i], tid[i]]
		
		self.data_header = [obj.sample for obj in self.quants]
			
		return

	def build_gene_id_level_tables(self):
		
		if self.num_quants < 1:
			error_message("No samples!")
			return None		
	
		# build union of gene id from each sample
		all_gid = set()
		for obj in self.quants:
			if obj.gid2idx is not None:
				all_gid.update(obj.gid2idx.keys())
		
		all_gid = sorted(list(all_gid))
		num_gid = len(all_gid)
		
		# now build the output
		tid = [None for i in range(num_gid)]		
		gnames = [None for i in range(num_gid)]
		
		self.counts = [[0 for i in range(self.num_quants)] for j in range(num_gid)]
		self.eff_counts = [[0 for i in range(self.num_quants)] for j in range(num_gid)]
		self.fpkm = [[0 for i in range(self.num_quants)] for j in range(num_gid)]
		self.tpm = [[0 for i in range(self.num_quants)] for j in range(num_gid)]
				
		for i in range(num_gid):
			gid = all_gid[i]
			for j in range(self.num_quants):
				obj = self.quants[j]
				if gid in obj.gid2idx:
					gidx = obj.gid2idx[gid]
						
					tid[i] = ";".join(obj.gid_annot[gid]['tid'])
					gnames[i] = ";".join(obj.gid_annot[gid]['gname'])
					
					self.counts[i][j] = obj.gid_hits[gidx]
					self.eff_counts[i][j] = obj.gid_eff_hits[gidx]
					self.fpkm[i][j] = obj.gid_fpkm[gidx]
					self.tpm[i][j] = obj.gid_tpm[gidx]
		
		self.annot = [[] for i in range(num_gid)]
		self.annot_header = ["gene_id", "gene_name", "transcript_id"]
		for i in range(num_gid):
			self.annot[i] = [all_gid[i], gnames[i], tid[i]]
		
		self.data_header = [obj.sample for obj in self.quants]
			
		return
		
	def build_report(self):
		
		# build initial table
		rtable = [[] for j in range(self.num_quants+1)]
		
		rtable[0] = ["Sample", "Reads", "Raw_MI", "MI", "Genes_Over_0_MI", "Genes_Over_1_MI"]
		if self.quants[0].has_sense:
			has_sense = True
		
		rtable[0] += ["q0", "q10", "q20", "q30", "q40", "q50", "q60", "q70", "q80", "q90"]
		
		for i in range(self.num_quants):
			obj = self.quants[i]
			j = i+1
			rtable[j] = [obj.sample, obj.total_reads, obj.total_raw_mi, obj.total_mi, obj.mi_over1, obj.mi_over2]
			for k in range(len(obj.bulk_bin_exp)):
				rtable[j].append(obj.bulk_bin_exp[k]*1.0/obj.mi_over1)
		
		return rtable

class Quant(object):
	
	def __init__(self):
		self.file = None
		self.sample = None
		
		#
		# relation tables
		#
		
		# transcript id to row index
		self.tid2idx = defaultdict(int)
		# gene id to transcript id
		self.gid2tid = defaultdict(set)
		# gene name to transcript id
		self.gn2tid = defaultdict(set)

		# raw rows from file
		self.annot = {}
		self.lengths = []
	
		self.hits = []
		self.eff_hits = []
		self.fpkm = []
		self.tpm = []		
		
		self.total_hits = 0
		self.total_eff_hits = 0
		self.total_fpkm = 0
		
		# for gene-level summaries
		self.gn2idx = defaultdict(int)
		self.gene_num_iso = None
		self.gene_exp_iso = None
		self.gene_annot = None
		self.gene_hits = None
		self.gene_eff_hits = None
		self.gene_fpkm = None
		self.gene_tpm = None

		# for gene_id-level summaries
		self.gid2idx = defaultdict(int)
		self.gid_num_iso = None
		self.gid_exp_iso = None
		self.gid_annot = None
		self.gid_hits = None
		self.gid_eff_hits = None
		self.gid_fpkm = None
		self.gid_tpm = None
	
	def parse(self, fname):
		self.file = fname
		
		tmp = fname.split("/")
		tmp2 = tmp[-1].split(".")
		self.sample = tmp2[0]
		
		# parse the file
		idx = 0
		with open(fname, "r") as fin:
			# read header to get column indices of annotation and lengths. after that
			# the columns should go hits, eff_hits, fpkm, tpm
			
			header = fin.readline().strip().split("\t")
			idx = 0
			found = False
			while idx < len(header):
				if header[idx]=="hits":
					found = True
					break
				self.annot[header[idx]] = []
				idx += 1
			
			if not found:
				error_message("Unable to locate the 'hits' column in {}".format(fname))
				sys.exit(1)
				
			annot_bound = idx
			
			# columns of data
			hit_col = annot_bound
			eff_hit_col = hit_col+1
			fpkm_col = eff_hit_col+1
			tpm_col = fpkm_col+1
			
			# parse on
			idx = 0
			for szl in fin:
				aln = szl.strip().split("\t")
				
				self.tid2idx[aln[0]] = idx
				idx += 1
				
				# copy annotation values
				for j in range(annot_bound):
					self.annot[header[j]].append(aln[j])
				
				self.gid2tid[aln[1]].add(aln[0])
				self.gn2tid[aln[2]].add(aln[0])
				
				self.hits.append(float(aln[hit_col]))
				self.eff_hits.append(float(aln[eff_hit_col]))
				self.fpkm.append(float(aln[fpkm_col]))
				self.tpm.append(float(aln[tpm_col]))
				
		
		self.total_hits = sum(self.hits)
		self.total_eff_hits = sum(self.eff_hits)
		self.total_fpkm = sum(self.fpkm)
		
		return

	def summarize_gene_level(self):
		
		# loop through gene names from gn2tid
		iso_count = 0
		iso_exp_count = 0
		
		self.gene_annot = {}

		num_genes = len(self.gn2tid.keys())

		self.gene_num_iso = [0 for i in range(num_genes)]
		self.gene_exp_iso = [0 for i in range(num_genes)]
		self.gene_hits = [0 for i in range(num_genes)]
		self.gene_eff_hits = [0 for i in range(num_genes)]
		self.gene_fpkm = [0 for i in range(num_genes)]
		self.gene_tpm = [0 for i in range(num_genes)]
		
		idx = 0
		for gn in self.gn2tid.keys():
			# hit up eacn isoform
			self.gn2idx[gn] = idx
			self.gene_annot[gn] = { 'tid':[], 'gid':set() }
			for tid in list(self.gn2tid[gn]):
				tid_idx = self.tid2idx[tid]

				#
				# deal with annotation of the gene
				#

				# concatenate tid, chrom, strand, length
				tid_hat = "{}|{}|{}|{}".format(tid, self.annot['chrom'][tid_idx], self.annot['strand'][tid_idx], 
					self.annot['length'][tid_idx])
				self.gene_annot[gn]['tid'].append(tid_hat)
				# concatenate gene id, chrom, strand
				gid_hat = "{}|{}|{}".format(self.annot['gene_id'][tid_idx], self.annot['chrom'][tid_idx], self.annot['strand'][tid_idx])
				self.gene_annot[gn]['gid'].add(gid_hat)
				
				#
				# collect data
				#
				self.gene_num_iso[idx] += 1
				if self.hits[tid_idx] > 5:
					self.gene_exp_iso[idx] += 1
				self.gene_hits[idx] += self.hits[tid_idx]
				self.gene_eff_hits[idx] += self.eff_hits[tid_idx]
				self.gene_fpkm[idx] += self.fpkm[tid_idx]
				self.gene_tpm[idx] += self.tpm[tid_idx]
				
				# done!
				
			# increment gene index
			idx += 1
		
		return


	def summarize_gene_id_level(self):
		
		# loop through gene names from gid2tid
		iso_count = 0
		iso_exp_count = 0
		
		self.gid_annot = {}

		num_gid = len(self.gid2tid.keys())

		self.gid_num_iso = [0 for i in range(num_gid)]
		self.gid_exp_iso = [0 for i in range(num_gid)]
		self.gid_hits = [0 for i in range(num_gid)]
		self.gid_eff_hits = [0 for i in range(num_gid)]
		self.gid_fpkm = [0 for i in range(num_gid)]
		self.gid_tpm = [0 for i in range(num_gid)]
		
		idx = 0
		for gid in self.gid2tid.keys():
			# hit up eacn isoform
			self.gid2idx[gid] = idx
			self.gid_annot[gid] = { 'tid':[], 'gname':set() }
			for tid in list(self.gid2tid[gid]):
				tid_idx = self.tid2idx[tid]

				#
				# deal with annotation of the gene id
				#

				# concatenate tid, chrom, strand, length
				tid_hat = "{}|{}|{}|{}".format(tid, self.annot['chrom'][tid_idx], self.annot['strand'][tid_idx], 
					self.annot['length'][tid_idx])
				self.gid_annot[gid]['tid'].append(tid_hat)
				# concatenate gene name, chrom, strand
				# gname_hat = "{}|{}|{}".format(self.annot['gene_name'][tid_idx], self.annot['chrom'][tid_idx], self.annot['strand'][tid_idx])
				gname_hat = self.annot['gene_name'][tid_idx]
				self.gid_annot[gid]['gname'].add(gname_hat)
				
				#
				# collect data
				#
				self.gid_num_iso[idx] += 1
				if self.hits[tid_idx] > 5:
					self.gid_exp_iso[idx] += 1
				self.gid_hits[idx] += self.hits[tid_idx]
				self.gid_eff_hits[idx] += self.eff_hits[tid_idx]
				self.gid_fpkm[idx] += self.fpkm[tid_idx]
				self.gid_tpm[idx] += self.tpm[tid_idx]
				
				# done!
				
			# increment gene index
			idx += 1
		
		return		

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Combine multiple expression quantifications from 'bam-gcounts-X' or 'bwa-quant'")
parser.add_argument('quants', type=str, nargs="+", 
	help="Specify 2 or more quantification files from 'bam-gcounts-X' for 'bwa-quant'")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

