#!/usr/bin/python
#==============================================================================
# salmon-to-isofrac.py
#
# Shawn Driscoll
# 20170512
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses Salmon based expression quantification and generates a per-gene id
# isoform fraction report.
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser

from collections import defaultdict

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")
HBIN = 16000

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	annot = gid2tid = gn2tid = None
	dquant = {}
	dktiso = {}
	
	#
	# figure out stub name from input file
	tmp = args.gtf.split("/")
	tmp2 = tmp[-1].split(".")
	stub = tmp2[0]
	
	
	#
	# load the GTF
	sys.stderr.write("loading gtf\n")
	annot, gid2tid, gn2tid = gtf_load(args.gtf)
	sys.stderr.write("done\n")

	#
	# load salmon file
	sys.stderr.write("Loading {}\n".format(args.ktiso))
	dktiso = load_ktiso(args.ktiso, annot, gid2tid)

	#
	# print report
	print "gene_id\ttranscript_id;fraction;length;eff_length"
	for gid in sorted(dktiso.keys()):
		lout = [gid]
		ltid = []
		for tid in sorted(dktiso[gid]['isoforms'].keys()):
			sz = "{};{};{};{}".format(tid, dktiso[gid]['isoforms'][tid]['fraction'], 
				dktiso[gid]['isoforms'][tid]['length'], dktiso[gid]['isoforms'][tid]['eff_length'])
			ltid.append(sz)

		lout.append("|".join(ltid))
		print "\t".join(lout)


		
	return 0


#==============================================================================
# general functions
#==============================================================================

def load_ktiso(fname, annot, gid2tid):
	# header 0:Name, 1:Length, 2:EffectiveLength, 3:TPM, 4:NumReads
	# header = { 'Name':0, 'Length':1, 'EffectiveLength':2, 'TPM':3, 'NumReads':4 }
	# we'll parse it in by gene id so we can bundle multi-isoform features. 
	# also we'll only pull in stuff from multi-isofrom features
	ddata = {}
	
	fin = open(fname, "r")
	szl = fin.readline()
	for szl in fin:
		aln = szl.strip().split("\t")
		tid = aln[0]
		
		if tid not in annot:
			warning("transcript {} is not in the GTF annotation".format(tid))
			continue
			
		hits = float(aln[4])
		gid = annot[tid]['gene_id']
			
		if gid not in ddata:
			ddata[gid] = {'hits':0, 'isoforms':{} }
		
		ddata[gid]['hits'] += hits
		ddata[gid]['isoforms'][tid] = { 'hits':hits, 'length':float(aln[1]), 'eff_length':float(aln[2]) }

	fin.close()
	
	# loop back through and calculate isoform fractions for multi-isoform genes
	for gid in ddata.keys():
		
		if len(gid2tid[gid])==1:
			# single isoform. give it 100% fraction to simplify 
			# downstream analysis
			for tid in ddata[gid]['isoforms'].keys():
				ddata[gid]['isoforms'][tid]['fraction'] = 1
		else:
			# get total hits
			t = ddata[gid]['hits']
			for tid in ddata[gid]['isoforms'].keys():
				t0 = ddata[gid]['isoforms'][tid]['hits']
				# calculate isoform fraction
				ddata[gid]['isoforms'][tid]['fraction'] = t0/t if t > 0 else 0
		
	return ddata	

def warning(sz):
	sys.stderr.write("Warning: {}\n")

#==============================================================================
# GTF functions
#==============================================================================

def gtf_parseline(sz):

	tmp = sz.strip().split("\t")

	grow = {
		'rname':tmp[0],
		'db':tmp[1],
		'type':tmp[2],
		'start':int(tmp[3]),
		'end':int(tmp[4]),
		'strand':tmp[6],
		'attrs':{}}

	# parse attributes
	fsplit = tmp[8].split("\"")
	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		grow['attrs'][key.strip()] = fsplit[i+1].strip()
		i += 2

	return grow

#
# print a parsed gtf row out in GTF format!
def gtfrow_as_gtf(grow):
	lout = [grow['rname'], 
		grow['db'], 
		str(grow['start']),
		str(grow['end']),
		".",
		grow['strand'], 
		"."
	]
	
	# build attribute string
	attr_string = "transcript_id \"{}\"; ".format(grow['attrs']['transcript_id'])
	attr_string += "gene_id \"{}\"; ".format(grow['attrs']['gene_id'])

	for atname in grow['attrs'].keys():
		if atname != "transcript_id" and atname != "gene_id":
			attr_string += "{} \"{}\";".format(atname, grow['attrs'][atname])
		
	lout.append(attr_string)
	
	return "\t".join(lout)
	

def gtf_transcript_id(grow):
	if "transcript_id" in grow['attrs']:
		return grow['attrs']['transcript_id']
	return None

def gtf_gene_id(grow):
	if "gene_id" in grow['attrs']:
		return grow['attrs']['gene_id']
	return None

def gtf_gene_name(grow):
	if "gene_name" in grow['attrs']:
		return grow['attrs']['gene_name']
	return None

def gtf_load(fname):
	#
	# load a gtf file into a few annotation tables

	dannot = {}
	dgid2tid = defaultdict(set)
	dgname2tid = defaultdict(set)

	fin = open(fname, "r")
	for szl in fin:
		grow = gtf_parseline(szl)

		tid = gtf_transcript_id(grow)
		if tid not in dannot:
			dannot[tid] = {
				'rname':grow['rname'],
				'strand':grow['strand'],
				'db':grow['db'],
				'gene_id':gtf_gene_id(grow),
				'gene_name':gtf_gene_name(grow),
				'num_exons':0,
				'start':grow['start'],
				'end':grow['end'],
				'length':0,
				'exons':[]}

		dannot[tid]['length'] += grow['end']-grow['start']+1
		dannot[tid]['exons'].append([grow['start'], grow['end']])
		dannot[tid]['num_exons'] += 1
		
		if grow['start'] < dannot[tid]['start']:
			dannot[tid]['start'] = grow['start']
		if grow['end'] > dannot[tid]['end']:
			dannot[tid]['end'] = grow['end']

		if dannot[tid]['gene_id'] is not None:
			dgid2tid[dannot[tid]['gene_id']].add(tid)

		if dannot[tid]['gene_name'] is not None:
			dgname2tid[dannot[tid]['gene_name']].add(tid)

	fin.close()

	#
	# pass through the annotation and sort the exons by position
	for tid in dannot.keys():
		dannot[tid]['exons'].sort(key=lambda x: x[0])

	return dannot, dgid2tid, dgname2tid

def gtf_num_exons(annot, tid):
	if tid not in annot:
		return -1
	
	return annot[tid]['num_exons']

def gtf_get_intron_chain(annot, tid):
	if tid not in annot:
		return None
	
	ints = []
	n = gtf_num_exons(annot, tid)
	if n < 2:
		return None

	elast = annot[tid]['exons'][0]
	for i in range(1, n):
		e = annot[tid]['exons'][i]
		ints.append([elast[1]+1, e[0]-1])
		elast = e
	
	return ints



def region_init(rname, start, end, tag=None, index=None):
	r = { 'rname':rname, 'start':start, 'end':end }
	
	if tag is not None:
		r['tag'] = tag

	if index is not None:
		r['index'] = index

	return r

#
# build a quick-lookup table for exon features from a loaded gtf
def gtf_lookup_table(annot):

	lktable = defaultdict(list)

	# loop through annotation by transcript id
	for tid in annot.keys():
		# loop through exons
		eidx = 0
		for e in annot[tid]['exons']:
			r = region_init(annot[tid]['rname'], e[0], e[1], tag=tid, index=eidx)
			h = region_hash(r)
			for hid in h:
				# insert the region in each binf
				lktable[hid].append(r)

	return lktable


#
# look up region r in the lookup table. return info from each hit including 
# the tag and index values
def gtf_find_hits(lktable, r):

	# hash the region, r
	lhash = region_hash(r)
	# hits list
	d_hits = {}

	# scan through hashs
	for h in lhash:
		if h in lktable:
			# scan through regions in this bucket
			for r0 in lktable[h]:
				rres = compare_regions(r, r0)
				if rres > 0:
					# we have a hit!
					ovl_len = region_overlap_length(r, r0)
					if r0['tag'] not in d_hits:
						d_hits[r0['tag']] = { 'index': [ r0['index'] ], 'length': [ovl_len[0]] }
					else:
						# add additional hit only if it is not to the same exon as one already
						# recorded. this would happen when we have a region that crosses 
						# the bin boundary and therefore would be found to hit the same feature
						# if that feature also crosses the bin boundary
						if r0['index'] not in set(d_hits[r0['tag']]['index']):
							d_hits[r0['tag']]['index'].append(r0['index'])
							d_hits[r0['tag']]['length'].append(ovl_len[0])

	#
	# return the dict of hits
	return d_hits



#
# assuming 'r' is a dict with 'rname', 'start' and 'end' fields
def region_hash(r):
	bin0 = binN = 0

	bin0 = int(r['start'])/HBIN
	binN = int(r['end'])/HBIN

	hout = ["{}:{}".format(r['rname'], bin0)]
	
	if binN > bin0:
		while bin0 < binN:
			bin0 += 1
			hout.append("{}:{}".format(r['rname'], bin0))

	return hout

#
# compare two regions. 
# return value:
# 0 for no overlap
# 1 for overlap
# 2 for identical
def compare_regions(r1, r2):

	# check ref names. if not equal then we're done
	if r1['rname'] != r2['rname']:
		return 0

	# ref names must be equal
	if r1['start']==r2['start'] and r1['end']==r2['end']:
		# starts and ends are identical
		return 2

	# now check for overlap
	if r1['end'] >= r2['start'] and r2['end'] >= r1['start']:
		# overlap!
		return 1

def region_length(r):
	return r['end']-r['start']+1

#
# calculate length of overlap between two regions.
# this is accomplished by finding the minimum value
# of 4 different measurements:
# A: length of r1
# B: length of r2
# C: end of r1 - start of r2
# D: end of r2 - start of r1
# if the minimum is negative then the result is 0
def region_overlap_length(r1, r2):
	len_A = region_length(r1)
	len_B = region_length(r2)
	len_C = r1['end']-r2['start']
	len_D = r2['end']-r1['start']

	rres = min([len_A, len_B, len_C, len_D])

	if rres <= 0:
		return 0

	# return length of overlap as well as ratios of the overlap to the length 
	# of the features
	return [ rres, rres*1.0/len_A, rres*1.0/len_B ]


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Convert a GTF to GFF3. Optionally only return genes with multiple isoforms.")
parser.add_argument('gtf', type=str, help="GTF related to the quantifications")
parser.add_argument('ktiso', type=str, 
	help="Isoform abundance estimates from Salmon")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


