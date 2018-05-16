#!/usr/bin/env python
#==============================================================================
# strip-un-and-mixed.py
# 
# Shawn Driscoll
# 20130408
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Strips out unaligned and mixed paired alignments from a SAM stream
#==============================================================================

import sys
import pysam

#
# variables
#

szl = ""
ll = []
argc = len(sys.argv)
argv = sys.argv[1:]

bamin = pysam.Samfile(argv[0], "rb")
bamout = pysam.Samfile("-", "wb", template=bamin)

#
# main script
#

for aln in bamin:

	if aln.flag & 0x4 or aln.flag & 0x8:
		bamout.write(aln)
	#}

#}


bamin.close()
bamout.close()
