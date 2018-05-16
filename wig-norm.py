#!/usr/bin/python
#==============================================================================
# wig-norm.py
#
# Shawn Driscoll
# 20170531
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Calculate sizes and normalize several wig files. We can optionally 
# also write out the binned versions of the wigs  
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
# used for printing exceptions
#import linecache
import traceback

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
HBIN = 16000

# region fields
REGION_RNAME = 0
REGION_START = 1
REGION_END = 2
REGION_TAG = 3

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	dbins = {}
	dlktable = {}
	lbins = []
	rhits = {}
	
	# calculate wig norm factors if more than one sample
	
	# this will be a matrix of the per region depths per file which 
	# we can then use to calculate normalization factors afterwards
	num_files = len(args.wig)
	depth_sets = [[] for i in range(num_files)]
	fidx = 0
	depths = [0 for i in range(num_files)]
	region_depths = [0 for i in range(num_files)]
	
	if args.size_table is None:
	
		#
		# build index
		#
		message("building genome index")
		dbins, lbins, dlktable = build_lookup_table(args)
		
		#
		# parse the file and assign hits
		#
		for fname in args.wig:
			foutname = fname + ".bin"
			
			message("intersecting {}".format(fname))	
			
			with open(fname, "r") as fin:
				clear_bins(dbins)
							
				for szl in fin:
					if re.search("^track", szl):
						continue
						
					aln = szl.strip().split("\t")
					r = region_init(aln[0], int(aln[1])+1, int(aln[2]))
					d = float(aln[3])
					depths[fidx] += d*region_length(r)
					
					# lookup hits
					h = region_hash(r)
					rhits = {}
					for hid in  h:
						if hid in dlktable:
							for r0 in dlktable[hid]:
								# overlap?
								rres = region_overlap_length(r, r0)
								if rres[0] > 0:
									# add in the depth of the region scaled by the size of the overlap
									# in bases
									#dbins[r0[REGION_TAG]] += d*rres[0]
									if r0[REGION_TAG] not in rhits:
										rhits[r0[REGION_TAG]] = d*rres[0]
					
					if len(rhits.keys()) > 0:
						for did in rhits.keys():
							dbins[did] += rhits[did]
							region_depths[fidx] += rhits[did]
							
			
			# write result to disk and populate this file's depth_set vector
			if args.write_bin:
				fout = open(foutname, "w")
				
			depth_sets[fidx] = [0 for i in range(len(lbins))]
			bin_idx = 0
			for did in lbins:
				if args.write_bin:
					lout = explode_region_string(did)
					lout[1] = str(int(lout[1])-1)
					lout.append("{}".format(dbins[did]))
					fout.write("\t".join(lout) + "\n")
					
				depth_sets[fidx][bin_idx] = dbins[did]
				bin_idx += 1
			
			fidx += 1
			
			if args.write_bin:
				fout.close()		
	
		if num_files > 1:
			# calculate normalization factors
			message("Calculating size factors for all samples")
			sizes = calc_norm_factors(depth_sets)
			# get median depth from region totals
			median_depth = nzmedian(region_depths)
			# depth per billion factor constat for all files
			dbp_factor = 1.0

			if args.dpb_norm:
				dbp_factor = 1000000000.0/median_depth
			
			# write size factors out along with depths
			
			with open("wig_norm_size_factors.tsv", "w") as fout:
				fout.write("sample\tsize_factor\tregion_depth\ttotal_depth\tdepth_ratio\n")
				for i in range(num_files):
					depth_ratio = region_depths[i]*100.0/depths[i]
					lout = [args.wig[i], "{}".format(sizes[i]), "{}".format(region_depths[i]), "{}".format(depths[i]), "{:0.2f}".format(depth_ratio)]
					fout.write("\t".join(lout) + "\n")
	
	else:
		# user provided a size table
		dsizes = defaultdict(float)
		region_depths = []
		dtmp = []
		with open(args.size_table, "r") as fin:
			for szl in fin:
				aln = szl.strip().split("\t")
				dsizes[aln[0]] = float(aln[1])
				dtmp.append(float(aln[2]))
		
		sizes = [-1 for i in range(num_files)]
		for i in range(num_files):
			if args.wig[i] in dsizes:
				sizes[i] = dsizes[args.wig[i]]
			else:
				warning_message("{} was not found in size table".format(args.wig[i]))
		
		# keep the region depths for those files that matched up to the size table
		for i in range(num_files):
			if sizes[i] > 0:
				region_depths.append(dtmp[i])
	
		# get median depth from region totals
		median_depth = nzmedian(region_depths)
		# depth per billion factor constat for all files
		dbp_factor = 1000000000.0/median_depth
	
	
	if num_files > 1:

		message("Creating normalized version of wigs")
		for i in range(num_files):
			
			if sizes[i] <= 0:
				warning_message("Skipping {} because it does not have a size factor".format(args.wig[i]))
				continue
			
			fname = args.wig[i]
			foutname = "{}_norm".format(fname)
			
			fin = open(fname, "r")
			fout = open(foutname, "w")
			for szl in fin:
				aln = szl.strip().split("\t")
				# scale depth to depth per billion 
				# and adjust by sample size factor
				d = float(aln[3])*dbp_factor/sizes[i]
				aln[3] = "{}".format(d)
				fout.write("\t".join(aln)+"\n")
			
			fout.close()
			fin.close()
		
		
	return 0

def mean(v):
	
	mu0 = 0
	mu = 0
	i = 0
	for x in v:
		i += 1
		mu0 = mu
		mu = mu0 + (x-mu0)/i
	
	return mu
	

def clear_bins(d):
	for k in d.keys():
		d[k] = 0

def explode_region_string(sz):
	tmp = sz.split(":")
	tmp2 = tmp[1].split("-")
	
	return [tmp[0]]+tmp2


def calc_norm_factors(x):
	
	xhat = transpose(x)
	n = len(xhat)
	xlog = []
	xmu = []
	
	# each "column" has to be log transformed. if there are zeros then 
	# we'll drop those now
	for i in range(n):
		# get mean and mean-centered log values
		mu, log_vals = geom_mean(xhat[i])
		if mu > 0:
			# keep it
			xlog.append(log_vals)
#			print mu, log_vals
	
	# now we transpose and find the medians. those are the size factors
	xloghat = transpose(xlog)
	m = len(xloghat)
	
	sizes = [1 for i in range(m)]
	for i in range(m):
		sizes[i] = math.exp(nzmedian(xloghat[i]))

	return sizes

def geom_mean(x):
	
	if min(x) <= 0:
		return 0, None
	
	xhat = []
	log_sum = 0
	for i in range(len(x)):
		xhat.append(math.log(x[i]))
	
	log_sum = sum(xhat)
	log_mu = log_sum/len(x)
	# subtract log mean from log values
	xhat_trans = [a-log_mu for a in xhat]
	
	return log_sum/len(x), xhat_trans
	

def transpose(x):
	m = len(x)
	n = len(x[0])
	
	xhat = []
	for i in range(n):
		xhat.append([0 for i in range(m)])
	
	for i in range(m):
		for j in range(n):
			xhat[j][i] = x[i][j]
	
	return xhat

def nzmedian(v):
	
	nnz = []
	n = len(v)
	med = None
	
	for i in range(n):
		if v[i] != 0:
			nnz.append(v[i])
		
	nhat = len(nnz)
	nmed = int(math.floor(nhat/2))
	
	tmp = sorted(nnz)
	nnz = tmp
	
	if (nhat % 2) != 0:
		# odd length
		med = (nnz[nmed]+nnz[nmed+1])*1.0/2
	else:
		med = nnz[nmed]
	
	return med

#
# build bins and lookup table. we also make a list of the bin keys
# so that we know what order to print them back out
def build_lookup_table(args):

	lktable = defaultdict(list)
	dbins = defaultdict(float)
	lbins = []
	
	if args.regions is not None:
		# parse bed file and build lookup table from that
		message("Loading regions from {}".format(args.regions))
		with open(args.regions, "r") as fin:
			for szl in fin:
				aln = szl.strip().split("\t")
				
				lbins.append(aln[3])
				dbins[aln[3]] = 0
				
				r = region_init(aln[0], aln[1], aln[2], tag=aln[3])
				h = region_hash(r)
				for hid in h:
					lktable[hid].append(r)
	
	else:	
		# open the chrom sizes file
		with open(args.chrom_sizes, "r") as fin:
			for szl in fin:
				aln = szl.strip().split("\t")
				rsize = int(aln[1])
				rname = aln[0]
				pos = 1
				while pos < rsize:
					right = pos+args.b-1
					if right > rsize:
						right = rsize
					
					# make a region
					rid = "{}:{}-{}".format(rname, pos, right)
					lbins.append(rid)
					dbins[rid] = 0
					
					r = region_init(rname, pos, right, tag=rid)
					h = region_hash(r)
					for hid in h:
						lktable[hid].append(r)
					
					pos += args.b

	return dbins, lbins, lktable


def region_init(rname, start, end, tag=None):
	# region is just a list
	r = [rname, int(start), int(end), None]
	# set optional fields
	if tag is not None:
		r[REGION_TAG] = tag
	
	return r

#
# assuming 'r' is a dict with 'rname', 'start' and 'end' fields
def region_hash(r):
	bin0 = binN = 0

	bin0 = int(r[REGION_START])/HBIN
	binN = int(r[REGION_END])/HBIN

	hout = ["{}:{}".format(r[REGION_RNAME], bin0)]
	
	if binN > bin0:
		while bin0 < binN:
			bin0 += 1
			hout.append("{}:{}".format(r[REGION_RNAME], bin0))

	return hout

#
# compare two regions. 
# return value:
# 0 for no overlap
# 1 for overlap
# 2 for identical
def compare_regions(r1, r2):

	rres = 0

	# check ref names. if not equal then we're done
	if r1[REGION_RNAME] != r2[REGION_RNAME]:
		return 0

	# ref names must be equal
	if r1[REGION_START]==r2[REGION_START] and r1[REGION_END]==r2[REGION_END]:
		# starts and ends are identical
		return 2

	# now check for overlap
	if r1[REGION_END] >= r2[REGION_START] and r2[REGION_END] >= r1[REGION_START]:
		# overlap!
		return 1
	
	return rres
	

def region_str(r):
	sz_pos = "{}:{}-{}".format(r[REGION_RNAME], r[REGION_START], r[REGION_END])

	if r[REGION_TAG] is not None:
		sz_pos += "|{}".format(r[REGION_TAG])


	return sz_pos

def region_length(r):
	return r[REGION_END]-r[REGION_START]+1

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
	len_C = r1[REGION_END]-r2[REGION_START]+1
	len_D = r2[REGION_END]-r1[REGION_START]+1

	rres = min([len_A, len_B, len_C, len_D])

	if rres <= 0:
		return [0, 0, 0]

	# return length of overlap as well as ratios of the overlap to the length 
	# of the features
	return [ rres, rres*1.0/len_A, rres*1.0/len_B ]


def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

#def PrintException():
#    exc_type, exc_obj, tb = sys.exc_info()
#    f = tb.tb_frame
#    lineno = tb.tb_lineno
#    filename = f.f_code.co_filename
#    linecache.checkcache(filename)
#    line = linecache.getline(filename, lineno, f.f_globals)
#    sys.stderr.write('EXCEPTION IN ({}, LINE {} "{}"): {}\n'.format(filename, lineno, line.strip(), exc_obj))
#    return

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


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('wig', nargs="+", metavar="wig", type=str, 
	help="Input wig or wig files to process")

parser.add_argument('-b', type=int, default=100000, help="Bin size in basepairs [100000]")
parser.add_argument('--write-bin', action="store_const", const=True, default=False, 
	help="Write the binned version of the wigs used for calculating size factors")
parser.add_argument('--no-norm', action="store_const", const=True, default=False, 
	help="Do not normalize the wig files. Use this in along with --write-bin to only export the intersection of data with regions")

parser.add_argument('--dpb-norm', action="store_const", const=True, default=False, 
	help="Apply 'depth-per-billion' scaling on top of the geometric mean based normalization. total depth is the median of all sample depths.")

parser.add_argument('--size-table', default=None, type=str, 
	help="Provide a table of size factors (as produced by this program). The first column must contain the file names.")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--regions', type=str, default=None, 
	help="BED format file of regions to use for normalization.")
group.add_argument('--chrom-sizes', type=str, default=None, 
	help="File containing the sizes of each chromosome. First column should be names and second should be sizes. Used if not using --regions")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
#	except KeyboardInterrupt:
#		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)
