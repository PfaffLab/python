#!/usr/bin/python
#
# numcols.py
#
# Shawn Driscoll
# 20140509
#
# pulls head of file or pipped lines from stdin and returns a column count 
# by tab delims per row
#

import sys,os
import subprocess as sp

try:

	argc = len(sys.argv)-1
	if argc == 0:
		sys.stderr.write("expecting input on stdin...\n")
		fin = sys.stdin
	else:
		p = sp.Popen("head {}".format(sys.argv[1]).split(), stdout=sp.PIPE)
		fin = p.stdout

	for szl in fin:
		aln = szl.strip().split("\t")
		print len(aln)

	fin.close()

except KeyboardInterrupt:
	sys.stderr.write("\n*cancelled*\n")


