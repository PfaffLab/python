#!/usr/bin/env python
#
# testing parsing sam alignments
#

import sys
import pysam as ps
#import pickle
from os.path import isfile
import re
from collections import defaultdict
from time import time
from GenomeJunk import PysamTools as ps_tools
from Basics import messages as ms

#==============================================================================
# other
#==============================================================================

	
class Stats(object):
	
	def __init__(self):
		self.lines = 0
		self.reads = 0
		self.primary_aligned = 0
		self.secondary_alignments = 0
		
		self.first_mate = 0
		self.second_mate = 0
		self.first_mate_secondary = 0
		self.second_mate_secondary = 0
		
		self.first_strand = 0
		self.second_strand = 0
		self.first_mate_first_strand = 0
		self.first_mate_second_strand = 0
		self.second_mate_first_strand = 0
		self.second_mate_second_strand = 0
	

def usage():
	sz = "usage: bam-stats.py <sam/bam>"
	return sz


#==============================================================================
# main
#==============================================================================

if __name__ == "__main__":

	argv = sys.argv
	argc = len(sys.argv)

	if argc < 2:
		ms.message(usage())
		sys.exit(1)

	argv = argv[1:len(argv)]

	stat_list = []

#	t0 = time()
#	lktable2 = ta.build_lktable()
#	ms.time_diff(t0)

	for i in range(len(argv)):
		if not isfile(argv[i]):
			ms.error_message("file does not exist ({})".format(argv[i]))
			continue
		
		mstat = Stats()		
	
		ms.message("parsing {}".format(argv[i]))
		t0 = time()
		with ps.AlignmentFile(argv[i]) as fin:
			rnames = ps_tools.get_alignmentfile_rnames(fin)
			numhit = 0
	
			for aln in fin:
				mstat.lines += 1
				
				if (aln.flag & 0x100) == 0:
					mstat.reads += 1
				
				if mstat.lines % 1000000 == 0:
					ms.progress_message("parsed {} lines".format(mstat.lines))
				
				if aln.flag & 0x4:
					continue
					
				if aln.flag & 0x100:
					mstat.secondary_alignments += 1
				else:
					mstat.primary_aligned += 1
				
				if aln.flag & 0x10:
					mstat.second_strand += 1
				else:
					mstat.first_strand += 1
				
				if aln.flag & 0x1:
					# paired end. collect additional info
					if aln.flag & 0x40:
						
						if aln.flag & 0x100:
							mstat.first_mate_secondary += 1
						else:
							mstat.first_mate += 1
						
						if aln.flag & 0x10:
							mstat.first_mate_second_strand += 1
						else:
							mstat.first_mate_first_strand += 1

					if aln.flag & 0x80:
						
						if aln.flag & 0x100:
							mstat.second_mate_secondary += 1
						else:
							mstat.second_mate += 1
						
						if aln.flag & 0x10:
							mstat.second_mate_second_strand += 1
						else:
							mstat.second_mate_first_strand += 1
			
			stat_list.append(mstat)
			ms.progress_message("parsed {} lines".format(mstat.lines), last=True)
			
		ms.time_diff(t0)
	
	print "\t".join(["file", "lines", "reads", "primary_aligned", "secondary_alignments", "first_strand", "second_strand", "first_mate", "second_mate", "first_mate_secondary", "second_mate_secondary", "fmfs", "fmss", "smfs", "smss"])
	for i in range(len(stat_list)):
		m = stat_list[i]
		lout = [argv[i], m.lines, m.reads, m.primary_aligned, m.secondary_alignments, 
			m.first_strand, m.second_strand, m.first_mate, m.second_mate, m.first_mate_secondary, m.second_mate_secondary, 
			m.first_mate_first_strand, m.first_mate_second_strand, m.second_mate_first_strand, m.second_mate_second_strand]
		
		print "\t".join(map(str, lout))
