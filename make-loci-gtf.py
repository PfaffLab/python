#!/usr/bin/env python
#==============================================================================
# make-loci-gtf.py
#
# shawn driscoll
# 20120809
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Using a transcript ID map (maybe created using gffread from cufflinks) this 
# script produces a loci summarized GTF
#==============================================================================

import sys,argparse
import pybedtools
import HTSeq

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Using a transcript ID map (maybe created using gffread from cufflinks) this script produces a loci summarized GTF")
parser.add_argument('gtf',type=str,help="GTF to compress into loci")
#parser.add_argument('map',type=str,help="Loci map where each line contains transcript IDs to merge into single loci")
args = parser.parse_args()

# check for file
try:
	fin = open(args.gtf)
except IOError as e:
	print "I/O error({0}): {1}".format(e.errno,e.strerror)
	sys.exit()
else:
	fin.close()


#==============================================================================
# parseGtf
# parse GTF annotation into a dict keyed by transcript ids. Each element contains
# both a list of the exon features for the transcript as well as a GenomicInterval
# object that spans the full range of the transcript.
#==============================================================================
def parseGtf(sz_file,d_g):

	##
	# variables
	szTid = ""

	# create reader for GTF
	gr = HTSeq.GFF_Reader(sz_file)

	# read through the GTF file and load into dict
	for feature in gr:
		if feature.type == "exon":
			szTid = feature.attr['transcript_id']
			feature.name = szTid

			# add to dict
			if szTid not in d_g:
				d_g[szTid] = {'exons':[], 'feature':None}

			# append current feature to transcript's exon list
			# kill strand
			feature.iv.strand = "."
			d_g[szTid]['exons'].append(feature)

			# update feature's genomic interval
			if d_g[szTid]['feature'] is None:
				d_g[szTid]['feature'] = HTSeq.GenomicFeature(feature.name,"gene",feature.iv.copy())
				d_g[szTid]['feature'].attr = feature.attr.copy()
			else:
				d_g[szTid]['feature'].iv.start = min(d_g[szTid]['feature'].iv.start,feature.iv.start)
				d_g[szTid]['feature'].iv.end = max(d_g[szTid]['feature'].iv.end,feature.iv.end)

#==============================================================================
# getParsedGtfTranscriptFeatures
# Returns a list of all of the GenmicFeature objects stored in d_g, a parsed
# GTF in dict from parseGtf (above)
#==============================================================================
def getParsedGtfTranscriptFeatures(d_g):

	# create empty list
	l_out = []

	# loop through features and load into list
	for szTid in d_g:
		l_out.append(d_g[szTid]['feature'])

	return l_out

#==============================================================================
# ivMergeSort
# Mergesort for GenomicFeature objects.
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
# Merge segment of mergesort function
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

			if ll[nil].iv.chrom < lr[nir].iv.chrom:
				l_result.append(ll[nil])
				nil += 1
			elif ll[nil].iv.chrom > lr[nir].iv.chrom:
				l_result.append(lr[nir])
				nir += 1
			else:
				# chrom values are equal so compare the start positions
				if ll[nil].iv.start <= lr[nir].iv.start:
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
# featureToBedLine
# Returns a tab delim string BED format representation of a GenomicFeature
# object.
#==============================================================================
def featureToBedLine(feature):
	return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(feature.iv.chrom,feature.iv.start,feature.iv.end,feature.name,".",feature.iv.strand)

#==============================================================================
# featureToGTFLine
# Returns a GTF tab delim string of a genomic feature
#==============================================================================
def featureToGtfLine(feature):
	
	l_out = [
		feature.iv.chrom,
		"loci",
		"exon",
		str(feature.iv.start+1),
		str(feature.iv.end),
		"0.000000",
		feature.iv.strand,
		".",
		"gene_id \"{0}\"; transcript_id \"{1}\";".format(feature.attr['gene_id'],feature.attr['transcript_id'])
		]

	return "\t".join(l_out)



#==============================================================================
# collapseSortedGF
# Collapse a list of sorted genomic features by merging overlapping genomic
# feature intervals from adjacent elements. 
#==============================================================================
def collapseSortedGF(ll):
	
	gf_new = None
	l_out = []
	i = 0

	# loop through list
	for feature in ll:

		if gf_new is None:
			#gf_new = feature
			gf_new = HTSeq.GenomicFeature(feature.name,feature.type,feature.iv.copy())
			gf_new.attr = feature.attr.copy()
		else:
			if feature.iv.overlaps(gf_new.iv):
				# features overlap, merge them
				gfMerge(gf_new,feature)
			else:
				# features don't overlap so append the completed "new" feature
				# to the output list
				l_out.append(gf_new)
				# start new "new" feature
				#gf_new = feature
				gf_new = HTSeq.GenomicFeature(feature.name,feature.type,feature.iv.copy())
				gf_new.attr = feature.attr.copy()

	# append final "new" feature
	l_out.append(gf_new)

	return l_out

#==============================================================================
# intersectGfAttr
# For a list of GenomicFeatures, create a consensus attr set for 'transcript_id'
# and 'gene_id' then update each feature to have the same values
#==============================================================================
def intersectGfAttr(l_gf):
	
	l_tid = []
	l_gid = []

	# loop through list to create consensus sets
	for feature in l_gf:
		l_tid += feature.attr['transcript_id'].split(",")
		l_gid += feature.attr['gene_id'].split(",")

	s_tid = set(l_tid)
	s_gid = set(l_gid)

	# loop through list again to update values
	for feature in l_gf:
		feature.attr['transcript_id'] = ",".join(list(s_tid))
		feature.attr['gene_id'] = ",".join(list(s_gid))




#==============================================================================
# gfMerge
# Merge two genomic features by modifying the first of the two genomic
# feature arguments (g1)
#==============================================================================
def gfMerge(g1,g2):

	# unique names and ids for merged gf
	s_names = set((g1.name).split(",") + (g2.name).split(","))
	s_tids = set((g1.attr['transcript_id']).split(",") + (g2.attr['transcript_id']).split(","))
	s_gids = set((g1.attr['gene_id']).split(",") + (g2.attr['gene_id']).split(","))

	# update fields

	g1.attr['transcript_id'] = ",".join(list(s_tids))
	g1.attr['gene_id'] = ",".join(list(s_gids))

	# alter iv
	g1.iv.start = min(g1.iv.start,g2.iv.start)
	g1.iv.end = max(g1.iv.end,g2.iv.end)
	# update name
	g1.name = ",".join(list(s_names))



#==============================================================================
# main
#==============================================================================
def main(args):
	
	## 
	# variables
	dGtf = {}
	i = 0
	l_features = None
	l_sortedFeatures = None
	l_loci = None
	l_lociFeatures = None

	##
	# parse gtf file
	sys.stderr.write("parsing gtf...\n")
	parseGtf(args.gtf,dGtf)

	## 
	# first we'll create a list of transcripts (full end to end ranges) sorted
	# by start position. We'll iterate through the sorted list collapsing features
	# down to loci which will result in a list of transcript id's belonging
	# to common loci

	sys.stderr.write("sorting list of genomic features...\n")

	# get list of genomic features
	l_features = getParsedGtfTranscriptFeatures(dGtf)
	# sort
	l_sortedFeatures = ivMergeSort(l_features)
	# collapse the features
	sys.stderr.write("collapsing list of genomic features...\n")
	l_loci = collapseSortedGF(l_sortedFeatures)

	# loop through the loci to build individual lists of genomic features
	# from each transcript in each locus
	for feature in l_loci:
		l_tlist = (feature.name).split(",")
		l_features = []
		for szTid in l_tlist:
			l_features += dGtf[szTid]['exons']

		# sort and collapse
		l_sortedFeatures = ivMergeSort(l_features)
		l_lociFeatures = collapseSortedGF(l_sortedFeatures)

		# print GTF records
		#print "locus: ",featureToGtfLine(feature)

		# update info for each feature so they all match since not all exons in a 
		# locus match between all features
		intersectGfAttr(l_lociFeatures)

		# print GTF records
		for gf in l_lociFeatures:
			print featureToGtfLine(gf)


#==============================================================================
# entry point
#==============================================================================
if __name__ == "__main__":
	main(args)
