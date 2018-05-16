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

def mergeSort(lv):
	# if list is length 1 then it's sorted - sweet
	if len(lv) == 1:
		return lv

	# else we need to split the list in two and pass both halves to this function
	n_mid = len(lv)/2
	l_left = [lv[i] for i in range(n_mid)]
	l_right = [lv[i] for i in range(n_mid,len(lv))]

	# recurse
	l_left = mergeSort(l_left)
	l_right = mergeSort(l_right)

	# merge lists
	return merge(l_left,l_right)

def merge(ll,lr):
	l_result = []

	n_ll = len(ll)
	n_lr = len(lr)
	nir = 0
	nil = 0

	# while counters are less than the lengths of left and right vectors
	while nir < n_lr or nil < n_ll:
		# if both counters are not at the end of their vectors, compare
		# next value from each
		if nil < n_ll and nir < n_lr:
			if ll[nil] <= lr[nir]:
				# append left value
				l_result.append(ll[nil])
				nil += 1
			else:
				# append right value
				l_result.append(lr[nir])
				nir += 1
		elif nil < n_ll:
			# left isn't done
			l_result.append(ll[nil])
			nil += 1
		elif nir < n_lr:
			# right isn't done
			l_result.append(lr[nir])
			nir += 1

	# finished merging
	return l_result


##
# test
l_unsorted = list(np.random.rand(int(sys.argv[1])))
n_ts = time.time()
l_sorted = mergeSort(l_unsorted)
n_te = time.time()

print "Total sort time: {:.2f} seconds".format(n_te-n_ts)

ar_diff = np.array(l_sorted[1:]) - np.array(l_sorted[0:len(l_sorted)-1])
print ar_diff.min()

#print l_unsorted
#print l_sorted

print len(l_unsorted),len(l_sorted)



