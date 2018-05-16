#!/usr/bin/env python
#==============================================================================
# gtf-select-multiiso.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
from collections import defaultdict
from igraph import Graph, WEAK
from time import localtime

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	annot = gid2tid = gn2tid = None
	num_multiiso = 0
	
	gidx = 0
		
	#
	# load the GTF
	message("loading {}".format(args.gtf))
	gtf = gtf_parse_into_refstrand_bins(args.gtf, args.unstranded)
	message("done")
	
	for ref in sorted(gtf.keys()):
		for strand in sorted(gtf[ref].keys()):
			foo = cluster_bin(gtf[ref][strand])
			fook = foo.keys()
			message("processing {}:{}; {} features".format(ref, strand, len(fook)))

			for idx in fook:
				# print back out all of the transcripts if the cluster has more than 1 transcript				
				if len(foo[idx]) > 1:
					num_multiiso += 1
					for tid in foo[idx]:
						for grow in gtf[ref][strand][tid]: 
							print gtfrow_as_gtf(grow)
		
		message("printed {} multi-isoform genes".format(num_multiiso))
	

	return 0

def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def cluster_bin(lbin):
	
	lverts = [] 
	ledges = []
	dverts = {}
	vcount = 0
	gid_groups = []
	
	n = len(lbin)
	
	# first we build this dict of transcript ids bundled by 3/5' ends of 
	# exons.
	dnodes = defaultdict(set)
	for tid in lbin.keys():
		n = len(lbin[tid])
		for i in range(n):
			v1 = "5:{}".format(lbin[tid][i]['start'])
			v2 = "3:{}".format(lbin[tid][i]['end'])
			dnodes[v1].add(tid)
			dnodes[v2].add(tid)
	
	
	# now we can loop back through this and build a graph with edges connecting
	# the transcript ids bundled in the nodes
	for vid in dnodes.keys():
		ltid = list(dnodes[vid])
		# first make sure the transcript ids are all in the vertex dict
		for tid in ltid:
			if tid not in dverts:
				dverts[tid] = vcount
				lverts.append(tid)
				vcount += 1
				
		# now make the edges
		for i in range(1, len(ltid)):
			e = [dverts[ltid[i-1]], dverts[ltid[i]]]
			ledges.append( [dverts[ltid[i-1]], dverts[ltid[i]]])
	
	# vcount has the total number of vertices. build graph and cluster
	g = Graph(directed=False)
	g.add_vertices(vcount)
	g.add_edges(ledges)
	clust = g.clusters(mode=WEAK)
	
	# the clust.membership vector is how we group the transcripts
	mem = clust.membership
	idx = range(len(mem))
	gid_groups = split(lverts, mem)
	return gid_groups
		
	
#
# like the R 'split' function. returns a dict
def split(subject, factor):
	# loop through factor and split into a dict
	dout = {}
	for i in range(len(factor)):
		fid = "{}".format(factor[i])
		if fid not in dout:
			dout[fid] = []
		dout[fid].append(subject[i])

	return dout
	
#
# defs
#

#
# parses a gtf into bins of chromosome and strand
# there is a third bin per chrom for any features that do not have strand
def gtf_parse_into_refstrand_bins(fname, unstranded):
	dout = {}
	
	fin = open(fname, "r")
	
	for szl in fin:
		grow = gtf_parseline(szl)
		
		rname = gtfrow_rname(grow)
		strand = gtfrow_strand(grow)
		tid = gtfrow_tid(grow)	
		
		if unstranded:
			strand = "ns"
			
		if rname not in dout:
			if unstranded:
				dout[rname] = { "ns": {} }
			else:
				dout[rname] = { "+":{}, "-":{}, "ns": {} }
		
		if strand=="+" or strand=="-":
			
			if tid not in dout[rname][strand]:
				dout[rname][strand][tid] = []
			dout[rname][strand][tid].append(grow)
			
		else:
			if tid not in dout[rname]['ns']:
				dout[rname]['ns'][tid] = []
			dout[rname]['ns'][tid].append(grow)
	
	fin.close()
	return dout
		
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

def gtfrow_rname(grow):
	return grow['rname']

def gtfrow_strand(grow):
	return grow['strand']

#
# print a parsed gtf row out in GTF format!
def gtfrow_as_gtf(grow):
	lout = [grow['rname'], 
		grow['db'], 
		"exon",
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
			attr_string += " {} \"{}\";".format(atname, grow['attrs'][atname])
		
	lout.append(attr_string)
	
	return "\t".join(lout)
	

def gtfrow_tid(grow):
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


parser = argparse.ArgumentParser(description="Search a GTF for redundant features.")
parser.add_argument('gtf', type=str, help="GTF to process")
parser.add_argument('-s', '--stub', type=str, default="GGID", 
	help="Stub for gene ids [GGID]")
parser.add_argument('-u', '--unstranded', action="store_const", const=True, default=False, 
	help="Ignore strand of isoforms")
args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


