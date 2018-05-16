#!/usr/bin/python


import sys,time,re
import HTSeq as hts


g_exons = hts.GenomicArrayOfSets("auto", stranded=False)

# start time
n_tStart = time.time()

gr = hts.GFF_Reader(sys.argv[1])
for feature in gr:
	if feature.type == "exon":
		sz_name = feature.attr['transcript_id'] + ";" + feature.attr['gene_id']
		if "gene_name" in feature.attr:
			sz_name += ";" + feature.attr['gene_name']

		sz_name += ";" + feature.iv.chrom

		g_exons[feature.iv] += sz_name

		# record total lengths of featurea in order to calculate RPKM later on
		if sz_name not in dLengths:
			dLengths[sz_name] = 0
			dHits[sz_name] = 0

		dLengths[sz_name] += feature.iv.end-feature.iv.start


