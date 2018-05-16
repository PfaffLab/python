#!/usr/bin/env python
#
# tcount-paired.py
#


import sys


if __name__ == "__main__":
	##
	## do this
	##

	tlist = {}

	fin = sys.stdin
	lline = []
	szl = ""
	paired = False
	flag = 0
	hit_factor = 1

	for szl in fin:
		lline = szl.strip().split("\t")

		flag = int(lline[1])
		tid = lline[2]

		paired = (flag & 0x1) != 0

		if (flag & 0x4) or (flag & 0x100):
			# skip unaligned and secondary
			pass
		else:

			if paired:
				if flag & 0x8:
					hit_factor = 1.0
				else:
					hit_factor = 0.5
			else:
				hit_factor = 1.0

			if tid not in tlist:
				tlist[tid] = 0

			tlist[tid] += hit_factor


	for tid in sorted(tlist.keys()):
		print "%s\t%0.2f" % (tid, tlist[tid])



