#!/usr/bin/python
#==============================================================================
# merge-sort.py
#
# Shawn Driscoll
# 20120809
#
# Basic merge sort
#==============================================================================

import numpy as np
import time
import sys

def selectionSort(lv):
	# variables
	ni = 0 
	nj = 0
	n_len = len(lv)
	n_min = None

	l_result = [lv[i] for i in range(n_len)]

	while ni < (n_len-1):
		# find min
		nj = ni
		n_min = None
		while nj < n_len:
			if n_min == None:
				n_min = nj
			else:
				if l_result[nj] < l_result[n_min]:
					n_min = nj

			nj += 1

		# swap values 
		n_tmp = l_result[ni]
		l_result[ni] = l_result[n_min]
		l_result[n_min] = n_tmp
		# current position holds the minimum value in the list including the currnet
		# position to the end so it's done. increment to next position
		ni += 1

	return l_result



##
# test
l_unsorted = list(np.random.rand(int(sys.argv[1])))
n_ts = time.time()
l_sorted = selectionSort(l_unsorted)
n_te = time.time()

print "Total sort time: {:.2f} seconds".format(n_te-n_ts)

ar_diff = np.array(l_sorted[1:]) - np.array(l_sorted[0:len(l_sorted)-1])
print ar_diff.min()

#print l_unsorted
#print l_sorted

print len(l_unsorted),len(l_sorted)



