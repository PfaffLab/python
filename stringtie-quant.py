#!/usr/bin/python
#==============================================================================
# stringtie-quant.py
#
# Shawn Driscoll
# 20170622
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Quantify isoform level abundance with Stringtie. 
#==============================================================================

import sys
import argparse
import math
import re
import traceback
from os.path import isfile, expanduser
from os import system
from collections import defaultdict
from time import localtime
import hashlib
import subprocess as sp

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
	
	for fname in args.bam:
		rres = run_stringtie(fname, args)

	return 0


def run_stringtie(fname, args):
	
	message("Processing {}".format(fname))
	
	stub = basename(fname)
	tmp_dir = hashlib.md5(fname).hexdigest()
	t_data = "{}/t_data.ctab".format(tmp_dir)
	dexpr = {}
	
	#
	# build command
	cmd_parts = ['stringtie',
		fname, 
		'-G', args.gtf, 
		'-x', '-B', '-b', tmp_dir, 
		'-e', '-p', args.p]
	
	if args.fr:
		cmd_parts.append('--fr')
	elif args.rf:
		cmd_parts.append('--rf')
	
	if args.u:
		cmd_parts.append('-u')
	
	cmd_parts = map(str, cmd_parts)
	
	cmd = " ".join(cmd_parts)
	
	sys.stderr.write("[system] {}\n".format(cmd))
	
	with open("/dev/null", "w") as dnull:
		p1 = sp.Popen(cmd_parts, stdout=dnull)
		p1.wait()
	
	# check for output
	if not isfile(t_data):
		error_message("Cannot find results! {}".format(t_data))
	else:
		# load this file up and produce a reformatted version
		
		dstie = {}
		message("Creating gene level output")
		with open(t_data, "r") as fin:
			# skip header
			szl = fin.readline()
			for szl in fin:
				srow = StringtieRow()
				srow.parse(szl)
				
				if srow.gene_id not in dstie:
					dstie[srow.gene_id] = Gene()
				
				dstie[srow.gene_id].add_isoform(srow)
			
		with open("{}.stie_gene".format(stub), "w") as fout:
			
			fout.write("gene_id\tgene_name\tlocus\tstrand\tnum_isoforms\ttranscript_list\tisofraction_hits\tisofraction_fpkm\n")
			for gid in sorted(dstie.keys()):
				dstie[gid].calc_isofractions()
				fout.write("\t".join(map(str, dstie[gid].tolist())) + "\n") 
				
		message("Creating isoform level output")
		system("mv {} {}.stie".format(t_data, stub))
		
	# kill the folder
	system("rm -r {}".format(tmp_dir))
	
	return 0

# split a file name up to base stub
def basename(fname):
	tmp = fname.split("/")
	tmp2 = tmp[-1].split(".")
	return tmp2[0]

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))
	return 0

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))
	return 0

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))
	return 0

def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0


class StringtieRow(object):
	
	def __init__(self):
		self.t_id = 0
		self.chr = ""
		self.strand = ""
		self.start = 0
		self.end = 0
		self.t_name = ""
		self.num_exons = 0
		self.length = 0
		self.gene_id = ""
		self.gene_name = ""
		self.cov = 0
		self.fpkm = 0
		
		self.hit_fraction = 0
		self.fpkm_fraction = 0
	
	def parse(self, sz):
		aln = sz.strip().split("\t")
		self.t_id = int(aln[0])
		self.chr = aln[1]
		self.strand = aln[2]
		self.start = int(aln[3])
		self.end = int(aln[4])
		self.t_name = aln[5]
		self.num_exons = int(aln[6])
		self.length = int(aln[7])
		self.gene_id = aln[8]
		self.gene_name = aln[9]
		self.cov = float(aln[10])
		self.fpkm = float(aln[11])
		return 0

class Gene(object):
	
	def __init__(self):
		self.name = ""
		self.id = ""
		self.isoforms = []
		self.num_isoforms = 0
		self.start = 0
		self.end = 0
		self.chrom = ""
		self.strand = ""
			
	# add isoform as a stringtie row
	def add_isoform(self, srow):
		self.isoforms.append(srow)
		
		if self.num_isoforms==0:

			self.start = srow.start
			self.end = srow.end
			self.chrom = srow.chr
			self.strand = srow.strand
			self.name = srow.gene_name
			self.id = srow.gene_id

		else:
			
			if srow.start < self.start:
				self.start = srow.start
			if srow.end > self.end:
				self.end = srow.end
		
		# increment isoform count
		self.num_isoforms += 1
		
		# sort isoforms on length in descending order
		self.isoforms.sort(key=lambda x:x.length, reverse=True)
		
	def get_locus_bounds(self):		
		# make a string that represents the boundary of this gene locus
		sz = "{}:{}-{}".format(self.chrom, self.start, self.end)
		return sz
	
	def calc_isofractions(self):
		# get total fpkm and total 'hit' counts
		total_fpkm = float(0)
		total_hit = float(0)
		
		for i in range(self.num_isoforms):
			total_fpkm += self.isoforms[i].fpkm
			total_hit += self.isoforms[i].fpkm*self.isoforms[i].length
		
		for i in range(self.num_isoforms):
			if total_fpkm > 0:
				self.isoforms[i].fpkm_fraction = self.isoforms[i].fpkm/total_fpkm
			if total_hit > 0:
				self.isoforms[i].hit_fraction = self.isoforms[i].fpkm*self.isoforms[i].length/total_hit
		
		return
	
	# create list format of what would be written out to a line
	def tolist(self):
		lout = [self.id, self.name, self.get_locus_bounds(), self.strand, self.num_isoforms]
		
		tid_list = ["" for i in range(self.num_isoforms)]
		fpkm_fraction = [0 for i in range(self.num_isoforms)]
		hit_fraction = [0 for i in range(self.num_isoforms)]
		
		for i in range(self.num_isoforms):
			tid_list[i] = self.isoforms[i].t_name
			fpkm_fraction[i] = "{:0.2f}".format(self.isoforms[i].fpkm_fraction)
			hit_fraction[i] = "{:0.2f}".format(self.isoforms[i].hit_fraction)
		
		lout.append(";".join(tid_list))
		lout.append(";".join(fpkm_fraction))
		lout.append(";".join(hit_fraction))
		
		return lout


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


parser = argparse.ArgumentParser(description="Quantify isoform level FPKM expression from aligned reads vs a GTF. Useful for establishing isoform fraction levels at genes.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
parser.add_argument('gtf', type=str, help="GTF annotation to quantify vs")
parser.add_argument('bam', type=str, nargs="+", help="Coordinate sorted BAM alignments to genome (1 or more files)")

strand_group = parser.add_mutually_exclusive_group(required=False)
strand_group.add_argument('--rf', action="store_const", const=True, default=False, 
	help="assume stranded library fr-firststrand")
strand_group.add_argument('--fr', action="store_const", const=True, default=False, 
	help="assume stranded library fr-secondstrand")

parser.add_argument('-p', type=int, default=1, 
	help="number of threads (CPUs) to use")
parser.add_argument('-u', action="store_const", const=True, default=False, 
	help="no multi-=mapping correction")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

