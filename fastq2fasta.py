#!/usr/bin/python
#
# fastq2fasta.py
#
# Shawn Driscoll
# 20140513
#
# changes fastq to fasta
#

import sys, os, re
import subprocess as sp

# variables
rname = ""
seq = ""
line = 0

try:
	argc = len(sys.argv)-1
	if argc == 0:
		sys.stderr.write("expecting input on stdin...\n")
		fin = sys.stdin
	else:
		# is the data gzipped?
		if re.search("\.gz$", sys.argv[1]):
			p = sp.Popen("gunzip -c {}".format(sys.argv[1]).split(), stdout=sp.PIPE)
			fin = p.stdout
		else:
			fin = open(sys.argv[1], "r")

	for szl in fin:
		line = line+1
		if line == 1:
			rname = szl.strip()
			# substitute the @ for a >
			rname = re.sub("^@", ">", rname)
		elif line == 2:
			seq = szl.strip()
			print rname + "\n" + seq
		elif line == 4:
			line = 0

	fin.close()

except KeyboardInterrupt:
	sys.stderr.write("\n*cancelled*\n")


