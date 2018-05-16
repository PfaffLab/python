#!/usr/bin/python
#==============================================================================
# modify-quant-ktiso.py
#
# Shawn Driscoll
# 20170512
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses a KTISO expression file (meaning Salmon, Sailfish or Kallisto type
# quantification) and applies the isofrom-fraction estimates per multi-isoform
# gene to the quantification from something like bam-gcounts or bwa-quant.
# single isoform gene levels are passed through without change. if a multi-isoform
# gene was not picked up in Salmon it will be zeroed out in the output.
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
from time import localtime
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

OUT_TID = 0
OUT_GID = 1
OUT_GNAME = 2
OUT_CHROM = 3
OUT_STRAND = 4
OUT_LENGTH = 5
OUT_EFF_LENGTH = 6
OUT_HITS = 7
OUT_EFF_HITS = 8
OUT_FPKM = 9
OUT_TPM = 10

#==============================================================================
# main
#==============================================================================

def main2(args):
	
	dgid = {}
	dgid_iso = {}
	dtid2gid = {}
	
	lheader = []
	
	
	quantF = QuantFile()
	
	# parse the quant file
	with open(args.aggene, "r") as fin:
		# skip header
		szl = fin.readline()
		quantF.add_header(szl)
		
		for szl in fin:
			
			q = QuantRow(szl)
			
			gid = q.gene_id
			tid = q.transcript_id
			hits = q.hits
			
			quantF.add_line(q)
						
			if gid not in dgid:
				dgid[gid] = {"count":0, "tid":{}}
			
			dgid[gid]['count'] += hits
			dgid[gid]['tid'][tid] = hits
			dtid2gid[tid] = gid
	
	# load isoform level stuff
	if args.stringtie is not None:
		
		with open(args.stringtie, "r") as fin:
			szl = fin.readline()
			for szl in fin:
				srow = StringtieRow()
				srow.parse(szl)
				
				hits = srow.length * srow.fpkm
				
				if srow.gene_id not in dgid_iso:
					dgid_iso[srow.gene_id] = {"count": 0, "tid": {}, "ifrac": {}}
				
				dgid_iso[srow.gene_id]['count'] += hits
				dgid_iso[srow.gene_id]['tid'][srow.t_name] = hits
				dgid_iso[srow.gene_id]['ifrac'][srow.t_name] = 0
	
	else:
		message("wut?")
		return 1
	
	
	#
	# calculate the isoform fraction values from the isoform level quant file
	#
	for gid in dgid_iso.keys():
		dropped = True
		# we're going to drop very small iso-fractions (< args.tolerance)
		while dropped:
			dropped = False
			
			hits = dgid_iso[gid]['count']
			min_ifrac = 1
			min_tid = ""
				
			for tid in dgid_iso[gid]['tid']:
				if hits > 0:
					tmp = dgid_iso[gid]['tid'][tid]
					ifrac = tmp/hits
					if (ifrac > 0) and (ifrac < min_ifrac):
						min_ifrac = ifrac
						min_tid = tid
					dgid_iso[gid]['ifrac'][tid] = ifrac
								
			if min_ifrac < args.tolerance:
				# subtract hits assigned to the min feature
				dgid_iso[gid]['count'] -= dgid_iso[gid]['tid'][min_tid]
				# zero hits to min feature
				dgid_iso[gid]['tid'][min_tid] = 0
				dgid_iso[gid]['ifrac'][min_tid] = 0
				dropped = True
	
	# now we should have the isoform fractions worked out
	for gid in dgid.keys():
		if gid not in dgid_iso:
			continue
		
		if dgid_iso[gid]['count'] > 0:
		
			tid0 = dgid[gid]['tid'].keys()
			all_there = True
			for tid in tid0:
				if tid not in dgid_iso[gid]['tid']:
					all_there = False
					message("Warning: {} is missing {} in the isoform quantification".format(gid, tid))
			
			if not all_there:
				continue
			
			# apply isoform fractions
			hits = dgid[gid]['count']
			for tid in tid0:
				dgid[gid]['tid'][tid] = hits*dgid_iso[gid]['ifrac'][tid] 
	
	print quantF.total_hits
	
	# now...
	# update the quant file
	for q in quantF.lines:
		# OK
		q.hits = dgid[q.gene_id]['tid'][q.transcript_id]

	quantF.update_total_hits()
	print quantF.total_hits
	
	quantF.update_eff_hits()
	quantF.update_fpkm()
	quantF.update_tpm()
	
	# print it out
	with open("{}iso".format(args.aggene), "w") as fout:
		fout.write(quantF.header)
		for l in quantF.lines:
			fout.write(l.__str__() + "\n")
	
	
	return 0

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
	fout_name = "{}iso".format(args.aggene)
	
	
	#
	# load the GTF
	message("loading {}".format(args.gtf))
	annot, gid2tid, gn2tid = gtf_load(args.gtf)
	message("done")

	#
	# load ktiso
	if args.salmon is not None:
		message("Loading {}".format(args.salmon))
		dktiso = load_ktiso(args.salmon, annot, gid2tid)
	elif args.stringtie is not None:
		message("Loading {}".format(args.stringtie))
		dktiso = load_stringtie(args.stringtie, annot, gid2tid)
	
	# load quant
	message("Loading {}".format(args.aggene))
	dquant = load_quant(args.aggene, annot)
	
	message("Checking for genes to import iso-fraction results into")
	gid_good = set()
	for gid in dquant.keys():
		if gid in dktiso:
			# confirm all of the transcripts are present
			keep = True
			if dktiso[gid]['hits'] < 1:
				keep = False
			else:
				for tid in dquant[gid]['isoforms'].keys():
					if tid not in dktiso[gid]['isoforms']:
						keep = False
			
			if keep:
				gid_good.add(gid)
	
	message("Importing iso-fraction results into {} genes".format(len(gid_good)))
	
	# loop through dquant and 'fix' the transcript abundances
	lres = []
	for gid in sorted(dquant.keys()):
		
		# get total hits from the quant file for this gene_id
		total_hits = dquant[gid]['hits']
		
		# loop through transcripts of this gene id
		for tid in dquant[gid]['isoforms'].keys():
			
			iso_fraction = 1
			# check for iso-fraction info from the ktiso file
			if gid in gid_good:
				iso_fraction = dktiso[gid]['isoforms'][tid]['fraction']
			
			# calculate corrected fraction
			corrected_hits = iso_fraction*total_hits
			
#			if corrected_hits < args.tolerance:
#				corrected_hits = 0
				
			# add output line
			ltmp = [
				tid, annot[tid]['gene_id'], annot[tid]['gene_name'], 
				annot[tid]['rname'], annot[tid]['strand'], 
				annot[tid]['length']]
			
			if dquant[gid]['isoforms'][tid]['eff_length'] is not None:
				ltmp.append(dquant[gid]['isoforms'][tid]['eff_length'])
			else:
				# use the salmon effective length estimate (bogus)
				ltmp.append(dktiso[gid]['isoforms'][tid]['eff_length'])
			
			eff_len = ltmp[OUT_EFF_LENGTH]
			
			# append hits
			ltmp.append(corrected_hits)
			# append eff_hits
			if eff_len > 0:
				ltmp.append( corrected_hits * annot[tid]['length']*1.0/eff_len )
			else:
				ltmp.append(0.0)
			
			lres.append(ltmp)
	
	message("Calculating FPKM and TPM")
		
	# now we can re-calculate fpkm and tpm
	total_hits = 0
	total_fpkm = 0

	# get sum total effective hits
	for l in lres:
		total_hits += l[OUT_HITS]
	
	# calc FPKM
	
	# NOTE: effective length is used to describe the effective transcript length
	# of the observed hit count. effective hits is used to describe the estimated
	# hits based on the full transcript length. so when calculating FPKM we should 
	# use either observed count and effective length or effective count and full 
	# length. 
	
	for l in lres:
		hits = float(l[OUT_HITS])
		length = l[OUT_EFF_LENGTH]
		f = 0
		if length > 0:
			f = hits*1e9/(length*total_hits)
		l.append(f)
		total_fpkm += f

	message("Printing results to {}".format(fout_name))
	
	fout = open(fout_name, "w")
	# write header
	fout.write("transcript_id\tgene_id\tgene_name\tchrom\tstrand\tlength\teff_length\thits\teff_hits\tfpkm\ttpm\n")
	
	# calc tpm and print output
	for l in lres:
		f = l[OUT_FPKM]
		t = f*1e6/total_fpkm
		l.append(t)
		fout.write("\t".join(map(str, l)))
		fout.write("\n")
	
	fout.close()
	
	return 0


#==============================================================================
# general functions
#==============================================================================


def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def load_ktiso(fname, annot, gid2tid):
	# header
	header = { 'Name':0, 'Length':1, 'EffectiveLength':2, 'TPM':3, 'NumReads':4 }
	# we'll parse it in by gene id so we can bundle multi-isoform features. 
	# also we'll only pull in stuff from multi-isofrom features
	ddata = {}
	
	with open(fname, "r") as fin:
		szl = fin.readline()
		for szl in fin:
			aln = szl.strip().split("\t")
			tid = aln[0]
			
			if tid not in annot:
				continue
				
			hits = float(aln[4])
			gid = annot[tid]['gene_id']
				
			if gid not in ddata:
				ddata[gid] = {'hits':0, 'isoforms':{} }
			
			ddata[gid]['hits'] += hits
			ddata[gid]['isoforms'][tid] = { 'hits':hits, 'eff_length':float(aln[2]) }
	
	# loop back through and calculate ratios
	for gid in ddata.keys():
		
		t = ddata[gid]['hits']
		for tid in ddata[gid]['isoforms'].keys():
			t0 = ddata[gid]['isoforms'][tid]['hits']
			# calculate isoform fraction
			ddata[gid]['isoforms'][tid]['fraction'] = t0/t if t > 0 else 0
		
	return ddata	

class QuantRow(object):
	def __init__(self, sz):
		aln = sz.strip().split("\t")
		
		self.transcript_id = aln[0]
		self.gene_id = aln[1]
		self.gene_name = aln[2]
		self.chrom = aln[3]
		self.strand = aln[4]
		self.length = float(aln[5])
		self.eff_length = float(aln[6])
		self.hits = float(aln[7])
		self.eff_hits = float(aln[8])
		self.fpkm = float(aln[9])
		self.tpm = float(aln[10])
		
	def __str__(self):
		return "\t".join(self.to_list())

	def to_list(self):
		lout = [self.transcript_id, self.gene_id, self.gene_name, self.chrom, self.strand, 
			str(self.length), "{:0.4f}".format(self.eff_length), "{:0.4f}".format(self.hits), 
			"{:0.4f}".format(self.eff_hits), "{:0.4f}".format(self.fpkm), "{:0.4f}".format(self.tpm)]
		return lout

class QuantFile(object):
	
	def __init__(self):
		self.header = []
		self.lines = []
		
		self.num_lines = 0
		self.total_hits = 0
		self.total_eff_hits = 0
		self.total_fpkm = 0

	def add_header(self, sz):
		self.header = sz
		return 0
	
	# append a quant object line to the file
	def add_line(self, q):
		self.lines.append(q)
		self.num_lines += 1
		
		self.total_hits += q.hits
		self.total_eff_hits += q.eff_hits
		self.total_fpkm += q.fpkm
		
		return 0
	
	def update_total_hits(self):
		self.total_hits = 0
		for i in range(self.num_lines):
			self.total_hits += self.lines[i].hits
		return 0
	
	def update_eff_hits(self):
		self.total_eff_hits = 0
		for i in range(self.num_lines):
			if self.lines[i].eff_length > 0:
				self.lines[i].eff_hits = self.lines[i].hits*float(self.lines[i].length)/self.lines[i].eff_length
				self.total_eff_hits += self.lines[i].eff_hits
			else:
				self.lines[i].eff_hits = 0
		return 0
		
	def update_fpkm(self):
		self.total_fpkm = 0
		for i in range(self.num_lines):
			f = self.lines[i].eff_hits*1e9/(float(self.lines[i].length)*self.total_eff_hits)
			self.total_fpkm += f
			self.lines[i].fpkm = f
		return 0
	
	def update_tpm(self):
		for i in range(self.num_lines):
			self.lines[i].tpm = self.lines[i].fpkm*1e6/self.total_fpkm
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

def load_stringtie(fname, annot, gid2tid):
	# header
	header = { 'Name':0, 'Length':1, 'EffectiveLength':2, 'TPM':3, 'NumReads':4 }
	# we'll parse it in by gene id so we can bundle multi-isoform features. 
	# also we'll only pull in stuff from multi-isofrom features
	ddata = {}
	
	with open(fname, "r") as fin:
		szl = fin.readline()
		for szl in fin:
			srow = StringtieRow()
			srow.parse(szl)
			
			tid = srow.t_name
			
			if tid not in annot:
				continue
				
			# scale the fpkm value by the length. this removes the transcript length
			# normalization making the set of values for the gene id proportional
			# to counts
			hits = srow.fpkm*srow.length
			gid = annot[tid]['gene_id']
				
			if gid not in ddata:
				ddata[gid] = {'hits':0, 'isoforms':{} }
			
			ddata[gid]['hits'] += hits
			ddata[gid]['isoforms'][tid] = { 'hits':hits, 'eff_length':srow.length, 'fraction': 0 }
	
	# loop back through and calculate ratios
	for gid in ddata.keys():
		
		t = ddata[gid]['hits']
		
		for tid in ddata[gid]['isoforms'].keys():
			t0 = ddata[gid]['isoforms'][tid]['hits']
			# calculate isoform fraction
			ddata[gid]['isoforms'][tid]['fraction'] = t0/t if t > 0 else 0
		
	return ddata	


def load_quant(fname, annot):
	header = None
	ddata = {}
	
	fin = open(fname, "r")
	header = header_to_dict(fin.readline())
	
	idx_tid = 0
	idx_hits = header['hits']
	idx_length = -1 if 'length' not in header else header['length']
	idx_eff_length = -1 if 'eff_length' not in header else header['eff_length']
		
	for szl in fin:
		aln = szl.strip().split("\t")
		
		tid = aln[idx_tid]
		if tid not in annot:
			continue
			
		gid = annot[tid]['gene_id']
		hits = float(aln[idx_hits])
		eff_length = length = None
		
		if idx_length > 0:
			length = float(aln[idx_length])
		if idx_eff_length > 0:
			eff_length = float(aln[idx_eff_length])
		
		if gid not in ddata:
			ddata[gid] = { 'hits':0, 'isoforms':{} }
		
		ddata[gid]['hits'] += hits
		ddata[gid]['isoforms'][tid] = { 'hits':hits, 'length':length, 'eff_length':eff_length }
	
	fin.close()
	return ddata

def header_to_dict(sz):
	aln = sz.strip().split("\t")
	d = defaultdict(int)
	i = 0
	for k in aln:
		d[k] = i
		i += 1
	
	return d

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


parser = argparse.ArgumentParser(description="Load the isoform level estimates from a quantification file and apply those to another quantification file", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
parser.add_argument('gtf', type=str, help="GTF related to the quantifications")
parser.add_argument('aggene', type=str, help="Gene level estimates from bam-gcounts or bwa-quant")

type_group = parser.add_mutually_exclusive_group(required=True)
type_group.add_argument('--salmon', type=str, default=None, 
	help="Isoform abundance estimates from Salmon")
type_group.add_argument('--stringtie', type=str, default=None,
	help="Isoform abundance estimates from Stringtie")

parser.add_argument('-t', '--tolerance', type=float, default=0.01, action="store", 
	help="Minimum isofraction threshold. Isoforms quantified to have lower than this value are zero'd")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main2(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


