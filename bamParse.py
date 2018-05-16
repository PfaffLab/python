#!/usr/bin/env python
#
# testing parsing sam alignments
#

import sys
import pysam as ps
import pickle
from os.path import isfile
import re
from collections import defaultdict
from time import time
from GenomeJunk import PysamTools as ps_tools
from GenomeJunk import *
from Basics import messages as ms

#==============================================================================
# other
#==============================================================================

def usage():
	sz = "usage: bamParse.py <refFlat> <sam/bam>"
	return sz

def organize_hits(hits):

	d = {}
	s = {}
	e = {}
	h = []
	for r in hits:
		tid = r.tag
		alen = r.attr['overlap_length']

		if tid not in d:
			d[tid] = 0
			s[tid] = r.strand
			e[tid] = set()

		d[tid] += alen
		e[tid].add(r.attr['exon_id'])

	for tid in d.keys():
		h.append([tid, d[tid], s[tid], ",".join(map(str, list(e[tid])))])

	return h


#==============================================================================
# main
#==============================================================================

if __name__ == "__main__":

	argv = sys.argv
	argc = len(sys.argv)

	if argc < 3:
		message(usage())
		sys.exit(1)

	argv = argv[1:len(argv)]

	ms.message("Loading refflat into object")
	ta = TranscriptomeAnnotation.TranscriptomeAnnotation()
	t0 = time()
	ta.load_refflat(argv[0])
	ms.time_diff(t0)

	ms.message("building lookup table from object")
	t0 = time()
	lktable2 = ta.build_lktable()
	ms.time_diff(t0)

	outfile = "{}.t.sam".format(argv[1])
	
	out_rnames = ta.names
	out_lengths = [ta.d[k].length for k in out_rnames]

	ms.message("parsing alignments and intersecting them with annotation")
	t0 = time()
	with ps.AlignmentFile(argv[1]) as fin, ps.AlignmentFile(outfile, "w", reference_names=out_rnames, reference_lengths=out_lengths) as fout:
		rnames = ps_tools.get_alignmentfile_rnames(fin)
		numhit = 0

		for aln in fin:
			if aln.flag & 0x4:
				continue

			# get alignment chunks
			aln_r = ps_tools.bam_aln_chunks(aln)
			hits = []
			for r in aln_r:
				rhat = GenericRegion.GenericRegion(rname=rnames[aln.rname], start=r[0], end=r[1], strand=ps_tools.bam_aln_strand_string(aln))
				chunk_hits = lktable2.lookup(rhat)
				if len(chunk_hits) > 0:
					hits += chunk_hits

			if len(hits) > 0:

				numhit += 1

				hitsHat = organize_hits(hits)
				
				for hit in hitsHat:

					ahat = ps.AlignedSegment()
					ahat.query_name = aln.query_name
					#ahat.query_sequence = aln.query_sequence
					#ahat.query_qualities = aln.query_qualities
					ahat.reference_id = ta.get_id(hit[0])
					ahat.cigar = aln.cigar
					ahat.mapping_quality = aln.mapping_quality
					ahat.reference_start = 0
					ahat.tags = aln.tags
					#ahat.tags.append(("OV", hit[1]))
					ahat.set_tag("OV", hit[1], value_type="i")
					ahat.set_tag("EX", hit[3], value_type="Z")

					# set flag for relative to transcript strand
					if aln.is_reverse:
						# alignment is on reverse strand
						if hit[2] == "-":
							# same
							ahat.flag = 0
						else:
							# reversed
							ahat.flag = 16
					else:
						# alignment is on forward strand
						if hit[2] == "-":
							# reversed
							ahat.flag = 16
						else:
							# same
							ahat.flag = 0

					fout.write(ahat)



	ms.time_diff(t0)
	ms.message("found {} hits".format(numhit))

