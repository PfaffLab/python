#!/usr/bin/python
#
# fastq2fasta.py
#
# Shawn Driscoll
# 20140513
#
# changes fasta to fastq. may not be suitable for very long sequences. 
# inserts default III qualities.
#

import sys, os, re
import subprocess as sp

def write_fastq(rname, seq):
	quals = "".join(["I" for i in range(len(seq))])

	print "@" + rname
	print seq
	print "+"
	print quals


# variables
rname = ""
seqnum = 0
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
		szl = szl.strip()

		if re.search("^>", szl):
			# new sequence. do we already have one collected?
			seqnum = seqnum + 1
			if len(seq) > 0:
				# print this bad boy out
				write_fastq("read" + str(seqnum) + "_" + rname, seq)

			rname = re.sub("^>", "", szl)
			seq = ""
		else:
			seq = seq + szl


	fin.close()

	# print out the final dude
	seqnum = seqnum+1
	if len(seq) > 0:
		write_fastq("read" + str(seqnum) + "_" + rname, seq)

except KeyboardInterrupt:
	sys.stderr.write("\n*cancelled*\n")


