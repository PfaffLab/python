#!/usr/bin/env python
#
# fastq-num-reads.py
#
# Read counter for fastq data
#

import os
import sys
import subprocess as sp
import re
from math import floor

def useage():
	sz = "useage: fastq-num-reads.py fastq {additional files...}"
	return sz

if __name__ == "__main__":

	argc = len(sys.argv) - 1
	if argc < 1:
		sys.stderr.write(useage() + "\n")
		sys.exit(1)

	argv = sys.argv[1:len(sys.argv)]
	argc = len(argv)
	nlines = 0
	fname = ""
	p1 = None

	for i in range(argc):
		nlines = 0
		fname = argv[i]
		# is the file gzipped?
		if re.search("\.gz$", fname):
			# open with subprocess
			p1 = sp.Popen("gunzip -c {}".format(fname).split(), stdout=sp.PIPE)
			if not p1:
				sys.stderr.write("Unable to open gz file with subprocess\n")
				sys.exit(1)
			fin = p1.stdout

		else:
			fin = open(fname, "r")

		for szl in fin:
			nlines += 1

		fin.close()

		nlines = int(floor(nlines / 4))

		print "{}\t{}".format(fname, nlines)



