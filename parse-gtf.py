#!/usr/bin/env python
#==============================================================================
# parse-gtf
#
# Shawn Driscoll
# 20120503
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses a GTF file turning the features column into multiple columns - useful
# for parsing a GTF with other scripts.
#==============================================================================

import sys
import argparse
import os.path
from GenomeJunk import GtfTools as gt
from Basics import messages as ms
from time import time
import re

#==============================================================================
# main
#==============================================================================
def main(args):
	
	# check input file
	if not os.path.isfile(args.gtf):
		ms.error_message("Input file is missing")
		return 1
	
	t0 = time()
	ms.message("Loading annotation")
	dtid, dtid2gname, all_attrs, tid_order = load(args.gtf)
	ms.time_diff(t0)
	
	tmp = all_attrs.difference(set(["transcript_id", "gene_id", "gene_name"]))
	all_attrs = sorted(list(tmp))
	
	# header
	sys.stdout.write("chrom\tdb\tfeature\tstart\tend\tscore\tstrand\tframe\ttranscript_id\tgene_id\tgene_name\t")
	sys.stdout.write("\t".join(all_attrs))
	sys.stdout.write("\n")
	
	for tid in tid_order:
		for gg in dtid[tid]:
			# print each row
			lout = gg.parts[0:8]
			lout.append(gg.transcript_id())
			lout.append(gg.gene_id())
			lout.append(gg.gene_name())
			
			for aid in all_attrs:
				if aid in gg.attr:
					lout.append(gg.attr[aid])
				else:
					lout.append("na")
			
			sys.stdout.write("\t".join(lout) + "\n")
	
	return 0
	
	
#==============================================================================
# functions
#==============================================================================

def load(fname):
	
	dtid = {}
	dtid2gname = {}
	tid_order = []
	dattrs = set(["transcript_id", "gene_id", "gene_name"])
	gene_names = set()
	n = 0
	
	with open(fname, "r") as fin:
		
		for szl in fin:
			if szl[0]=="#":
				# comment line
				continue
			
			gg = gt.GtfRow(szl)
			
			if gg.feature() != "exon":
				continue
			
			n += 1
			if (n % 100000) == 0:
				ms.progress_message("parsed {} lines. {} transcripts in {} genes".format(n, len(dtid.keys()), len(gene_names)))
			
			tid = gg.transcript_id()
			
			dattrs.update(gg.attr.keys())
			
			if tid not in dtid:
				dtid[tid] = []
				tid_order.append(tid)
			
			dtid[tid].append(gg)
			
			if tid not in dtid2gname:
				dtid2gname[tid] = gg.gene_name()
				gene_names.add(gg.gene_name())
	
	ms.progress_message("parsed {} lines. {} transcripts in {} genes".format(n, len(dtid.keys()), len(gene_names)), last=True)
	
		
	
	return dtid, dtid2gname, dattrs, tid_order
	
	

#==============================================================================
# check arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Parse a GTF file into tabbed columns", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("gtf", action="store", type=str, help="GTF file or - to read from stdin")
parser.add_argument('-o', type=str, action="store", default=None, help="Output file name stub. Base name without extension is used otherwise.")
args = parser.parse_args()

#==============================================================================
# variables
#==============================================================================

#==============================================================================
# main script
#==============================================================================

if __name__ == "__main__":
	sys.exit(main(args))
	
	

