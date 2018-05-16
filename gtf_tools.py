#!/usr/bin/env python

import sys

_CHROM = 0
_START = 1
_END = 2
_STRAND = 3
_TYPE = 4
_ATTR = 5

#==============================================================================
# gtf_to_dict
#
# Parse a GTF file into a dict keyed by sz_key.
#==============================================================================
def gtf_to_dict(sz_file,sz_key="transcript_id"):

	#
	# variables
	fin = None
	szl = ""
	d_out = {}
	ll = []
	lattr = []
	dattr = []
	lnew = []
	i = 0

	#
	# try opening file
	try:
		fin = open(sz_file,"r")
	except IOError,e:
		sys.stderr.write("I/O error({0}): {1}\n".format(e.errno,e.strerror))
		sys.exit()

	#
	# continue
	for szl in fin:
		ll = szl.strip().split("\t")

		# parse attributes
		lattr = ll[8].split(";")
		dattr = {}
		i = 0
		while i < len(lattr)-1:
			ltmp = lattr[i].split("\"")
			dattr[ltmp[0].strip()] = ltmp[1].strip()
			i += 1

		# create single exon object
		lnew = []
		lnew.append(ll[0]) # chrom
		lnew.append(int(ll[3])-1) # start
		lnew.append(int(ll[4])) # end
		lnew.append(ll[6]) # strand
		lnew.append(ll[2]) # type
		lnew.append(dattr.copy()) # attributes dict
		
		# get id
		szId = dattr[sz_key]
		if szId not in d_out:
			# 
			# append exon to the appropriate id
			#
			d_out[szId] = []
		
		d_out[szId].append(lnew)

	fin.close()

	return d_out
		

#==============================================================================
# gtf_to_list
#
# Parse a GTF file into a list with one entry per feature
#==============================================================================
def gtf_to_list(sz_file,sz_type="exon"):

	#
	# variables
	fin = None
	szl = ""
	l_out = []
	ll = []
	lattr = []
	dattr = []
	lnew = []
	i = 0

	#
	# try opening file
	try:
		fin = open(sz_file,"r")
	except IOError,e:
		sys.stderr.write("I/O error({0}): {1}\n".format(e.errno,e.strerror))
		sys.exit()

	#
	# continue
	for szl in fin:
		szl = szl.strip()
		ll = szl.split("\t")

		if ll[2] == sz_type:

			# parse attributes
			lattr = ll[8].split(";")
			dattr = {}
			i = 0
			while i < len(lattr)-1:
				ltmp = lattr[i].split("\"")
				dattr[ltmp[0].strip()] = ltmp[1].strip()
				i += 1

			lnew = []
			lnew.append(ll[0]) # chrom
			lnew.append(int(ll[3])-1) # start
			lnew.append(int(ll[4])) # end
			lnew.append(ll[2]) # type
			lnew.append(ll[6]) # strand
			lnew.append(dattr.copy()) # attributes dict
			
			l_out.append(lnew)

	fin.close()
	
	return l_out

		
#
# index_sorted_gtf_list
# Build a dict index of a sorted (by chrom and start position) list of 
# GTF features. Index holds start and end indices for each chromosome of 
# features in lgtf.
#
def index_sorted_gtf_list(lgtf):

	d_index = {}
	n_index = 0
	szlast = ""
	for lf in lgtf:
		if lf[_CHROM] not in d_index:
			d_index[lf[_CHROM]] = []
			# start of current chromosome in list
			d_index[lf[_CHROM]].append(n_index)
			# mark the end of the last chromosome as well
			if szlast != "":
				d_index[szlast].append(n_index)

		szlast = lf[_CHROM]
		n_index += 1

	d_index[szlast].append(n_index)

	return d_index


#
# iv_merge_sort
# Mergesort a bunch of parsed GTF features in a list
#
def iv_merge_sort(lv):
	# if list is length 1 then it's sorted - sweet
	if len(lv) == 1:
		return lv

	# else we need to split the list in two and pass both halves to this function
	n_mid = len(lv)/2
	l_left = [lv[i] for i in range(n_mid)]
	l_right = [lv[i] for i in range(n_mid,len(lv))]

	# recurse
	l_left = iv_merge_sort(l_left)
	l_right = iv_merge_sort(l_right)

	# merge lists
	return iv_merge(l_left,l_right)

#
# iv_merge
# Merge segment of mergesort function
#
def iv_merge(ll,lr):
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

			if ll[nil][_CHROM] < lr[nir][_CHROM]:
				l_result.append(ll[nil])
				nil += 1
			elif ll[nil][_CHROM] > lr[nir][_CHROM]:
				l_result.append(lr[nir])
				nir += 1
			else:
				# chrom values are equal so compare the start positions
				if ll[nil][_START] <= lr[nir][_START]:
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



def copy_feature(lf):

	l_new = [lf[i] for i in range(len(lf)-1)]
	l_new.append(lf[_ATTR].copy())
	return l_new

def feature_overlap(lf1,lf2):

	if lf1[_START] == lf2[_START] and lf1[_END] == lf2[_END]:
		# 
		# features are equal
		#
		return 2
	
	if lf1[_START] < lf2[_END] and lf1[_END] > lf2[_START]:
		#
		# features overlap
		#
		return 1

	# nothing
	return 0


#
# mergeList
# Merge overlapping adjacent features to produce a list of uniqe non overlapping
# features.
#
def merge_list(ll):

	l_out = []
	l_new = None
	n_index = 0

	for lf in ll:

		if l_new is None:
			l_new = copy_feature(lf)
			l_new.append(n_index)
		else:
			# check for overlap
			if l_new[_START] < lf[_END] and l_new[_END] > lf[_START]:
				# found overlap

				# modify l_new's boundaries to encompass both features
				l_new[_START] = min(l_new[_START],lf[_START])
				l_new[_END] = max(l_new[_END],lf[_END])

				# modify l_new's attributes to include both
				for szKey in l_new[_ATTR]:
					s_new = set(l_new[_ATTR][szKey].split(","))
					s_new.update(set(lf[_ATTR][szKey].split(",")))
					l_new[_ATTR][szKey] = ",".join(list(s_new))
			else:
				# no overlap, append to output and start a new feature
				l_out.append(l_new)
				l_new = copy_feature(lf)
				l_new.append(n_index)

		n_index += 1

	# append final feature
	l_out.append(l_new)

	return l_out

	