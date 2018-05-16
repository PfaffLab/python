#!/usr/bin/env python
#==============================================================================
# sam-read.py
#
# Shawn Driscoll
# 20120815
#
# 
#==============================================================================

import sys,time,re


def main(argv):
	if len(argv) == 0:
		return 1

	n_ts = 0
	n_te = 0

	# open reads
	fin = open(argv[0],"r")

	# iterate through file
	lkeep = ("M","D","N")
	n_ts = time.time()
	for szl in fin:
		szl = szl.strip()
		ll = szl.split("\t")
		# iterate through cigar

		l_features = []

		sz_strand = "+"
		if int(ll[1]) & 0x10:
			sz_strand = "-"

		# parse cigar line
		l_op = re.split("[0-9]+",ll[5])[1:]
		l_len = map(int,re.split("[A-Z]",ll[5])[0:len(l_op)])

		opKeep = ["M","D","N"]

		n_start = int(ll[3])
		n_end = n_start

		for i in range(len(l_op)):
			if l_op[i] in opKeep:
				n_end = n_start + l_len[i]

			if l_op[i] == "M":
				l_features.append([ll[2],n_start,n_end,sz_strand])

			n_start = n_end

		#print l_features

	fin.close()

	n_te = time.time()
	print "total time: {:.2e} seconds".format(n_te-n_ts)

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))

