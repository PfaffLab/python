#!/usr/bin/env python
#==============================================================================
# nexpr.py
#
# Shawn Driscoll
# 20120813
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script accepts a BAM file and  GTF annotation file to compute read counts
# and FPKM values for transcripts and loci specified in the GTF file. I've 
# put my own genomic feature pipeline together that seems to outperform 
# HTSeq's tools by nearly 400%. The GTF file is parsed into a hash where each
# hash key is the chromosome (or reference) name concatenated with each
# feature's start position divided by some bin value.
#==============================================================================

import sys,argparse,re
import math
#import HTSeq
import pysam
import time

# the following provides a slight performance increase
import gc
gc.disable()

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Quantify hits and expression in RPKM/FPKM from a SAM file against some reference annotation (GTF format)")
parser.add_argument("alignments",type=str,help="Alignments in SAM format or - to read from STDIN")
parser.add_argument("gtf",type=str,help="Gene annotation in GTF format")
parser.add_argument("-q",type=int,dest="n_qmin",default=0,help="Minimum mapping quality (MAPQ) from alignments. (default: 0)")
parser.add_argument("-v",type=int,dest="n_maxmm",default=-1,help="Maximum mismatches accepted for quantification. (default: any)")
parser.add_argument("--allow-singletons",dest="b_allowSingles",action="store_const",const=True,default=False,help="For paired reads. Allow reads that were not properly paired to be quantified. (default: False)")
parser.add_argument("--total-hits-norm",dest="b_totalHitsNorm",action="store_const",const=True,default=False,help="Normalize FPKM by all accepted alignments. (default: transcriptome alignments only)")
parser.add_argument("--strict",dest="b_strict",action="store_const",const=True,default=False,help="Strict mode. Requires 0x100, 0x200 and 0x400 to not be set in SAM lines before processing further")
parser.add_argument("--mf-out",dest="sz_mfout",action="store",type=str,default="",help="Write reads that hit multiple features out to this file.")
parser.add_argument("--stranded",dest="b_stranded",action="store_const",const=True,default=False)
parser.add_argument("-o",dest="sz_prefix",type=str,default="samexpr",help="Prefix for output files. (default: samexpr)")
parser.add_argument("--no-loci",dest="b_noLoci",action="store_const",const=True,default=False,help="Skip loci detection/quantification. Useful for annotations based on EST libraries or already compressed loci annotations.")

args = parser.parse_args()

#==============================================================================
# globals
#==============================================================================

_CHROM = 0
_START = 1
_END = 2
_STRAND = 3
_TYPE = 4
_ATTR = 5

HASH_BIN = 12000

#==============================================================================
# main
# 
#==============================================================================
def main(args):

	##
	# variables
	n_ts = 0
	n_te = 0
	l_gtf = []
	l_gtfSorted = []
	dStats = {'accepted':0, 'hits':0, 'skipped':0, 'multi-aligned':0, 'multi-feature':0}
	l_comp = []
	d_index = {}
	fin = None
	sf = None
	n_hits = 0
	n_multi = 0
	n_read = 0
	b_accepted = False
	b_paired = False
	n_hitfactor = 1
	l_loci = []

	#n_ts = time.time()
	#l_gtf = gtfToList(args.gtf)

	# parse in gtf data as list and as a dict
	n_ts = time.time()
	l_gtf,d_gtf = readGtf(args.gtf)
	
	if not args.b_noLoci:
		# need to create loci level info so we have to extract a list of the overall feature
		# boundaries and merge them down into unique regions
		l_featureBoundaries = []
		for szTid in d_gtf.keys():
			l_featureBoundaries.append(copyFeature(d_gtf[szTid]['feature']))

		# sort
		l_sorted = ivMergeSort(l_featureBoundaries)
		# compress
		l_loci = mergeList(l_sorted)

	n_te = time.time()
	sys.stderr.write("Loaded {:d} features ({:d} transcripts in {:d} loci) in {:.2f} seconds\n".format(len(l_gtf),len(d_gtf),len(l_loci),n_te-n_ts))

	#n_te = time.time()
	#sys.stderr.write("GTF parsed in {:.4f} seconds\n".format(n_te-n_ts))

	#n_ts = time.time()
	#l_gtfSorted = ivMergeSort(l_gtf)
	#n_te = time.time()
	#sys.stderr.write("Sort time: {:.4f} seconds\n".format(n_te-n_ts))

	d_hash = hashGtfList(l_gtf)

	#l_comp = mergeList(l_gtfSorted)

	# index
	#d_index = indexSortedGTFList(l_comp)

	#
	# open alignments file and find overlaps
	try:
		fin = open(args.alignments,"r")
	except IOError,e:
		sys.stderr.write("I/O error({0}): {1}\n".format(e.errno,e.strerror))
		return e.errno
	else:
		fin.close()

	# create pysam object to read from bam file
	sf = pysam.Samfile(args.alignments,"rb")
	fin = sf.fetch()

	# start timer
	n_ts = time.time()

	n_accepted = 0

	# continue
	for aln in fin:
		#szl = szl.strip()

		n_read += 1

		#b_accepted = alignmentAccepted(szl,args.b_strict,args.b_allowSingles,args.n_qmin,args.n_maxmm)
		b_accepted = alignmentAccepted(aln,args.b_strict,args.b_allowSingles,args.n_qmin,args.n_maxmm)

		# update paired status
		b_paired = aln.flag & 0x1

		if n_read % 2**17 == 0:
			n_te = time.time()
			sys.stderr.write("\rprocessed {:d} alignments ({:.2f} alignments/second)".format(n_read,float(n_read)/(n_te-n_ts)))

		if b_accepted:
			n_accepted += 1
			dStats['accepted'] += 1

			# get a list of matching genomic ranges from the sam alignment
			l_features = []
			l_features = parseSamLine(aln,sf)

			#print l_features
			# for each feature returned find overlaps
			s_hits = None
			n_tf = time.time()
			for lf in l_features:
				#lh = findOverlap(lf,l_comp,l_gtfSorted,d_index)
				lh = hashLookup(lf,d_hash)
				#fakeOverlap(lf,l_gtfSorted,d_index)
			
				if s_hits is None:
					s_hits = lh
				else:
					s_hits.update(lh)

			#if n_read % 100 == 0:
			#	print n_read,len(l_features),"{:.8f}".format(time.time()-n_tf),l_features

			n_hitFactor = 1
			if len(s_hits) > 0:
				n_hitFactor = 1.0/len(s_hits)

				for szId in list(s_hits):
					d_gtf[szId]['hits'] += n_hitFactor
				
				dStats['hits'] += 1
				#print szl
				#print s_hits

			if len(s_hits) > 1:
				dStats['multi-feature'] += 1

		else:
			dStats['skipped'] += 1

	
	n_te = time.time()
	sys.stderr.write("\rprocessed {:d} alignments ({:.2f} alignments/second\n".format(n_read,float(n_read)/(n_te-n_ts)))

	#for szId in sorted(d_gtf.keys()):
	#	print szId,d_gtf[szId]['hits']

	n_total = dStats['accepted'] + dStats['skipped']

	sys.stderr.write("alignments are ")
	if b_paired:
		sys.stderr.write("paired-end (only first read of pair is assigned)")
	else:
		sys.stderr.write("single-end")
	sys.stderr.write("\n")

	sys.stderr.write("{0} total alignments parsed; of these:\n".format(n_total))

	sys.stderr.write("  {:d} accepted ({:.2%})\n".format(dStats['accepted'],dStats['accepted']/float(n_total)))
	sys.stderr.write("  {:d} assigned ({:.2%})\n".format(dStats['hits'],dStats['hits']/float(n_total)))
	sys.stderr.write("  {:d} duplicated reads ({:.2%})\n".format(dStats['multi-aligned'],dStats['multi-aligned']/float(n_total)))
	sys.stderr.write("  {:d} assigned to multiple features ({:.2%})\n".format(dStats['multi-feature'],dStats['multi-feature']/float(n_total)))


	# print expression results to stdout
	if args.b_totalHitsNorm:
		n_norm = dStats['accepted']/1000000.
	else:
		n_norm = dStats['hits']/1000000.

	fout = open(args.sz_prefix + ".transcripts.expr","w")

	for szId in sorted(d_gtf.keys()):
		# calculate expression
		n_counts = round(d_gtf[szId]['hits'])
		n_expr = n_counts/(d_gtf[szId]['length']/1000.)/n_norm

		# print to output file
		l_out = [
			d_gtf[szId]['feature'][_ATTR]['gene_id'],
			d_gtf[szId]['feature'][_ATTR]['transcript_id'],
			"{0}:{1}-{2}".format(d_gtf[szId]['feature'][_CHROM],d_gtf[szId]['feature'][_START],d_gtf[szId]['feature'][_END]),
			"{:d}".format(d_gtf[szId]['length']),
			"{:d}".format(int(round(n_counts))),
			"{:.4f}".format(n_expr)]

		fout.write("\t".join(l_out) + "\n")
		#fout.write("\t".join(szId.split(";")) + "\t" + str(dLengths[szId]) + "\t" + str(int(math.floor(dHits[szId]))) + "\t" + str(n_expr) + "\n")

	fout.close()		

	if not args.b_noLoci:

		# prepare loci level counts and expression
		fout = open(args.sz_prefix + ".genes.expr","w")
		l_out = []
		for feature in l_loci:
			l_tids = feature[_ATTR]['transcript_id'].split(",")
			# sum counts from each tid
			n_counts = 0
			n_expr = 0
			for szTid in l_tids:
				n_counts += round(d_gtf[szTid]['hits'])
				#n_expr += dHits[dTranslation[szTid]]/(dLengths[dTranslation[szTid]]/1000.)/n_norm
				n_expr += round(d_gtf[szTid]['hits'])/(d_gtf[szTid]['length']/1000.)/n_norm

			# print to output file
			l_out = [
				feature[_ATTR]['gene_id'],
				feature[_ATTR]['transcript_id'],
				"{0}:{1}-{2}".format(feature[_CHROM],feature[_START],feature[_END]),
				"{:d}".format(int(round(n_counts))),
				"{:.4f}".format(n_expr)]

			fout.write("\t".join(l_out) + "\n")
			#fout.write("{0}\t{1}\n".format(feature.attr['gene_id'],feature.attr['transcript_id']))

		fout.close()


# 13512249 7405393

#==============================================================================
# readGtf
#
# Parse a GTF to prepare it for this analysis. This prepares a list of 
# every GTF feature parsed into a list for containing key information as well
# as a dict keyed by transcript_id. Each slot in the dict is a dict with a slot
# for hit counts, feature length in bases and a genomic feature to hold
# the overall boundary of the transcript. 
#==============================================================================
def readGtf(sz_file,sz_type="exon",szId="transcript_id"):

	#
	# variables
	fin = None
	szl = ""
	l_out = []
	ll = []
	lattr = []
	dattr = []
	lnew = []
	i = 0

	# features dict used to hold a unique entry for each feature keyed
	# by szId. this will be used to store hit counts and lengths
	d_features = {}

	#
	# try opening file
	try:
		fin = open(sz_file,"r")
	except IOError,e:
		sys.stderr.write("I/O error({0}): {1}\n".format(e.errno,e.strerror))
		sys.exit()

	#
	# continue
	for szl in fin:
		szl = szl.strip()
		ll = szl.split("\t")

		if ll[2] == sz_type:

			# parse attributes - don't know why i didn't think of this earlier
			lattr = ll[8].split(";")
			# pass attributes to a dict
			dattr = {}
			i = 0
			while i < len(lattr)-1:
				ltmp = lattr[i].split("\"")
				dattr[ltmp[0].strip()] = ltmp[1].strip()
				i += 1

			# create new list of key information
			lnew = []
			lnew.append(ll[0]) # chrom
			lnew.append(int(ll[3])-1) # start
			lnew.append(int(ll[4])) # end
			lnew.append(ll[2]) # type
			lnew.append(ll[6]) # strand
			lnew.append(dattr.copy()) # attributes dict

			# append feature to the list
			l_out.append(lnew)

			# prepare dict entry for this transcript if it doesn't exist already
			if dattr[szId] not in d_features:
				d_features[dattr[szId]] = {"hits":0,"length":0,"feature":None}

			# update length
			d_features[dattr[szId]]['length'] += lnew[_END]-lnew[_START]

			# update feature overall boundary
			if d_features[dattr[szId]]['feature'] is None:
				d_features[dattr[szId]]['feature'] = copyFeature(lnew)
			else:
				# update 
				d_features[dattr[szId]]['feature'][_START] = min(d_features[dattr[szId]]['feature'][_START],lnew[_START])
				d_features[dattr[szId]]['feature'][_END] = max(d_features[dattr[szId]]['feature'][_END],lnew[_END])

	# close input file
	fin.close()
	
	# return list and dict
	return (l_out,d_features)


#==============================================================================
# ivMergeSort
# Mergesort for genomic features stored in a list. 
#==============================================================================
def ivMergeSort(lv):
	# if list is length 1 then it's sorted - sweet
	if len(lv) == 1:
		return lv

	# else we need to split the list in two and pass both halves to this function
	n_mid = len(lv)/2
	l_left = [lv[i] for i in range(n_mid)]
	l_right = [lv[i] for i in range(n_mid,len(lv))]

	# recurse
	l_left = ivMergeSort(l_left)
	l_right = ivMergeSort(l_right)

	# merge lists
	return ivMerge(l_left,l_right)

#==============================================================================
# ivMerge
# Merge segment of mergesort function for genomic features. Features are 
# sorted by chrom and then start position.
#==============================================================================
def ivMerge(ll,lr):
	l_result = []

	n_ll = len(ll)
	n_lr = len(lr)
	nir = 0
	nil = 0

	# while counters are less than the lengths of left and right vectors
	while nir < n_lr or nil < n_ll:
		# if both counters are not at the end of their vectors, compare
		# next value from each
		if nil < n_ll and nir < n_lr:

			if ll[nil][_CHROM] < lr[nir][_CHROM]:
				l_result.append(ll[nil])
				nil += 1
			elif ll[nil][_CHROM] > lr[nir][_CHROM]:
				l_result.append(lr[nir])
				nir += 1
			else:
				# chrom values are equal so compare the start positions
				if ll[nil][_START] <= lr[nir][_START]:
					# append left value
					l_result.append(ll[nil])
					nil += 1
				else:
					# append right value
					l_result.append(lr[nir])
					nir += 1
		elif nil < n_ll:
			# left isn't done
			l_result.append(ll[nil])
			nil += 1
		elif nir < n_lr:
			# right isn't done
			l_result.append(lr[nir])
			nir += 1

	# finished merging
	return l_result


#==============================================================================
# featureOverlap
# Returns true if lf1 overlaps lf2. Simple test. If lf1's start is less than
# lf2's end and lf1's end is greater than lf2's start then the features overlap.
#==============================================================================
def featureOverlap(lf1,lf2):

	if lf1[_START] < lf2[_END] and lf1[_END] > lf2[_START]:
		return True

	return False


#==============================================================================
# mergeList
# Merge overlapping adjacent features to produce a list of unique non overlapping
# features. Each unique feature contains an extra index to the non-collapsed
# list for where each collapsed region starts.
#==============================================================================
def mergeList(ll):

	l_out = []
	l_new = None
	n_index = 0

	for lf in ll:

		if l_new is None:
			l_new = copyFeature(lf)
			l_new.append(n_index)
		else:
			# check for overlap
			if l_new[_START] < lf[_END] and l_new[_END] > lf[_START]:
				# found overlap

				# modify l_new's boundaries to encompass both features
				l_new[_START] = min(l_new[_START],lf[_START])
				l_new[_END] = max(l_new[_END],lf[_END])

				# modify l_new's attributes to include both
				for szKey in l_new[_ATTR]:
					s_new = set(l_new[_ATTR][szKey].split(","))
					s_new.update(set(lf[_ATTR][szKey].split(",")))
					l_new[_ATTR][szKey] = ",".join(list(s_new))
			else:
				# no overlap, append to output and start a new feature
				l_out.append(l_new)
				l_new = copyFeature(lf)
				l_new.append(n_index)

		n_index += 1

	# append final feature
	l_out.append(l_new)

	return l_out

#==============================================================================
# copyFeature
# Return a copy of the list feature, lf. Necessary since Python creates 
# references of list objects when assigning with =
#==============================================================================
def copyFeature(lf):

	l_new = [lf[i] for i in range(len(lf)-1)]
	l_new.append(lf[_ATTR].copy())
	return l_new


#==============================================================================
# alignmentAccepted
# Returns True if an alignment should be considered for analysis based on 
# filtering options specified at the command line.
#==============================================================================
def alignmentAccepted(aln,b_strict,b_allowSingles,n_qmin,n_maxmm):
	#
	# variables
	n_qual = aln.mapq
	n_flag = aln.flag
	n_mm = 0
	
	if (b_strict and ((n_flag & 0x700) == 0)) or not b_strict:
		# process further

		# check quality
		if n_qual >= n_qmin:
			# find reported mismatches

			if n_maxmm >= 0:
				n_mm = 0
				# get mismatches from options
				try:
					n_mm = aln.opt("XM")
				except KeyError,e:
					n_mm = aln.opt("NM")

			if n_maxmm < 0 or n_mm <= n_maxmm:

				# process paired or single end
				if n_flag & 0x1:
					# paired-end alignments
					if b_allowSingles:
						if n_flag & 0x40:
							return True

					else:
						if n_flag & 0x40 and n_flag & 0x2:
							return True

				else:
					# single-end alignments
					return True

	return False


#==============================================================================
# parseSamLine
# Parses cigar information and returns a set of features for any "M" features.
#==============================================================================
def parseSamLine(aln,sf):
	#
	# variables
	l_features = []

	sz_strand = "+"
	if aln.flag & 0x10:
		sz_strand = "-"

	# cigar ops to keep (see SAM format spec)
	ckeep = (0,2,3)

	# start a alignment start position (1 based)
	n_start = aln.pos+1
	n_end = n_start

	# loop through cigar string
	for ctype,clen in aln.cigar:
		# move end position only when M, D or N are encountered
		if ctype in ckeep:
			n_end = n_start + clen

		# if type is M then append feature to return list
		if ctype == 0:
			l_features.append([sf.getrname(aln.rname),n_start,n_end,sz_strand])

		# move start to end value (sometimes these are already equal)
		n_start = n_end

	# return list of features
	return l_features


#==============================================================================
# getFeatureHash
# Hash function for hashing the chromosome and position index of features.
#==============================================================================
def getFeatureHash(chrom,val):
	return chrom + str(val / HASH_BIN)

#==============================================================================
# hashGtfList
# Creates a dict from the list of GTF features in l_features. Each feature
# is hashed (with getFeatureHash) and that hash is used as a key in the 
# dict. The result is a dict with 1000's of bins, each containing several 
# GTF features. This allows super fast lookup when searching for overlaps since
# we can jump straight to a bin by using the same hash function on a query 
# feature.
# with the mm9 known gene GTF this yeilds bins with 5 to 15 features (mean of 8).
# so on average a lookup will have only  8 iterations.
#==============================================================================
def hashGtfList(l_features):

	#
	# variables
	dHash = {}
	lf = None
	n_sbin = 0
	n_ebin = 0

	# loop through features
	for lf in l_features:
		# find bin indicies. I added this workaround when I discovered an exon
		# that's nearly 100kb long in the mouse genome. It would span multiple 
		# bins (at the default 16kb bin length).
		n_sbin = lf[_START]/HASH_BIN
		n_ebin = lf[_END]/HASH_BIN

		# feature will span more than one bin so it will be added to EACH
		# in full.
		if n_sbin != n_ebin:

			# loop across bins
			i = n_sbin
			while i <= n_ebin:
				# manually assemble "hash"
				shash = lf[_CHROM] + str(i)
				# insert key if it doesn't exist
				if shash not in dHash:
					dHash[shash] = []

				# append feature to hash key
				dHash[shash].append(copyFeature(lf))

				i += 1

		else:
			# feature does not span multiple bins, assign to a single bin
			shash = lf[_CHROM] + str(n_sbin)

			if shash not in dHash:
				dHash[shash] = []

			dHash[shash].append(copyFeature(lf))

	# sort each hash position
	#for szh in dHash.keys():
	#	l_sorted = ivMergeSort(dHash[szh])
	#	dHash[szh] = [l for l in l_sorted]

	return dHash

#==============================================================================
# hashLookup
# Search for feature overlaps of lf with features stored in the search hash,
# dhash. By "hashing" lf we can jump straight to a small bin, or pair of bins,
# that would contain any possible overlaps.
#==============================================================================
def hashLookup(lf,dhash):

	# crete hash for feature's start and end
	shash = getFeatureHash(lf[_CHROM],lf[_START])
	ehash = getFeatureHash(lf[_CHROM],lf[_END])

	# if neither hash exists then this alignment doesn't have any overlap
	if shash not in dhash and ehash not in dhash:
		return set([])

	# start empty list for hits
	l_hits = []

	# if both the start and end hashes are identical then this feature's
	# potential overlaps are in one bin. if not then overlaps may be
	# in two bins. 
	if shash == ehash:
		if shash in dhash:
			l_features = dhash[shash]
			n_features = len(l_features)
			i = 0
			#while i < n_features and l_features[i][_START] < lf[_END]:
			while i < n_features:	
				if featureOverlap(l_features[i],lf):
					l_hits.append(l_features[i][_ATTR]['transcript_id'])

				i += 1

			#sys.stderr.write("{} iterations\n".format(i))
	else:
		# search both bins. more than two is unlikely unless the alignment length
		# surpasses the size of the bins (16kb)
		if shash in dhash:
			l_features = dhash[shash]
			n_features = len(l_features)
			i = 0
			#while i < n_features and l_features[i][_START] < lf[_END]:
			while i < n_features:	
				if featureOverlap(l_features[i],lf):
					l_hits.append(l_features[i][_ATTR]['transcript_id'])

				i += 1

			#sys.stderr.write("{} iterations\n".format(i))
		
		if ehash in dhash:
			l_features = dhash[ehash]
			n_features = len(l_features)
			i = 0
			#while i < n_features and l_features[i][_START] < lf[_END]:
			while i < n_features:	
				if featureOverlap(l_features[i],lf):
					l_hits.append(l_features[i][_ATTR]['transcript_id'])

				i += 1

			#sys.stderr.write("{} iterations\n".format(i))

	return set(l_hits)	



#==============================================================================
# main entry point
#==============================================================================
if __name__ == "__main__":
	sys.exit(main(args))
