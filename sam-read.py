#!/usr/bin/env python
#==============================================================================
# sam-read.py
#
# Shawn Driscoll
# 20120815
#
# 
#==============================================================================

import sys,pysam,time


def main(argv):
	if len(argv) == 0:
		return 1

	n_ts = 0
	n_te = 0

	# open bam file
	sf = pysam.Samfile(argv[0],"rb")

	# get iterator
	isf = sf.fetch()

	# iterate through file
	ckeep = (0,2,3)
	n_ts = time.time()
	for aln in isf:
		# iterate through cigar
		
		try:
			n_mm = aln.opt("XM")
		except KeyError,e:
			n_mm = aln.opt("NM")

		print n_mm
		
		l_features = []

		sz_strand = "+"
		if aln.flag & 0x10:
			sz_strand = "-"

		n_start = aln.pos+1

		for ctype,clen in aln.cigar:
			if ctype in ckeep:
				n_end = n_start + clen

			if ctype == 0:
				l_features.append([sf.getrname(aln.rname),n_start,n_end,sz_strand])

			n_start = n_end

		#print l_features

	n_te = time.time()
	print "total time: {:.2e} seconds".format(n_te-n_ts)

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))

