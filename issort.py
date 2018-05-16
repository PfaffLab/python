#!/usr/bin/env python
#==============================================================================
# issort.py
#
# Shawn Driscoll
# 20120826
#
# Insertion sort
#==============================================================================

import sys,random

def issort(data, compare):
	
	i = 0
	j = 0
	tmp = None
	data_size = len(data)
	
	for j in range(1,data_size):
		# copy data at j to tmp
		tmp = data[j]
		
		i = j - 1
		
		# determine position at which to insert the element. iterate backwards
		# as long as the data at i is greater than the data at tmp
		while i >= 0 and compare(data[i],tmp) > 0:
			# move data to the right to make space for current tmp data
			data[i+1] = data[i]
			i -= 1
		
		# i is now at the insertion position. it may be -1 since that would be 
		# allowed by the previous loop
		data[i+1] = tmp
		
	return 0

def qksort_partition(data, i, k, compare):
	
	pval = None
	temp = None
	
	# set of values for median of three
	lr = [(random.random() % (k - i + 1)) + i for n in range(3)]
	issort(lr,compare_int)
	
	# pivot index is center index in lr
	pval = data[int(lr[1])]
	
	#
	# create two partitions around pval
	i -= 1
	k += 1
	
	while True:
		while True:
			k -= 1
			if compare(data[k],pval) <= 0:
				break
		
		while True:
			i += 1
			if compare(data[i],pval) >= 0:
				break
		
		if i > k:
			break
		else:
			# swap data
			temp = data[i]
			data[i] = data[k]
			data[k] = temp
	
	return k

def qksort(data, i, k, compare):
	j = 0
	result = 0
	
	while i < k:
		j = qksort_partition(data,i,k,compare)
		if j < 0:
			# bail out
			return -1
		
		# recursivly sort left partition
		result = qksort(data,i,j,compare)
		if result < 0:
			# bail out
			return -1
		
		# iterate and sort right
		i = j + 1
	
	return 0

def compare_int(i,j):
	if i < j:
		return -1
	
	if i > j:
		return 1
	
	return 0
	
def main():
	
	lnums = [random.random() for i in range(20)]
	
	print lnums
	
	#issort(lnums,compare_int)
	qksort(lnums,0,len(lnums)-1,compare_int)
	
	print lnums
	

if __name__ == "__main__":
	sys.exit(main())