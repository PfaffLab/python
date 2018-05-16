#!/usr/bin/env python

import sys,time
import numpy as np

#
# Test HTSeq's GTF parsing speed against doing it manually

#gc.disable()

try:
	fin = open(sys.argv[1])
except IOError,e:
	print e
else:
	fin.close()

_CHROM = 0
_START = 1
_END = 2
_STRAND = 3
_TYPE = 4
_ATTR = 5


def parseGtfLine(szl):
	ll = szl.strip().split("\t")

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
	lnew.append(ll[6]) # strand
	lnew.append(ll[2]) # type
	lnew.append(dattr.copy()) # attributes dict	

	# return
	return lnew

def getFeatureHash(chrom,val,n_bin=16000):
	return chrom + str(val/n_bin)

#==============================================================================
# hashGtf
#
# Parse a GTF into a hash to make for quick lookups based on genomic
# location
#==============================================================================
def hashGtf(sz_file,n_bin=16000):

	# 
	# variables
	fin = None
	szl = ""
	dHash = {}

	try:
		fin = open(sz_file,"r")
	except IOError,e:
		sys.stderr.write("I/O error({0}): {1}".format(e.errno,e.strerror))
		sys.exit()

	# continue
	for szl in fin:
		lg = parseGtfLine(szl)

		sz_hash = getFeatureHash(lg[_CHROM],lg[_START],n_bin=n_bin)

		if sz_hash not in dHash:
			dHash[sz_hash] = []

		dHash[sz_hash].append(copyFeature(lg))

	fin.close()

	# sort each bucket
	for sz_hash in dHash.keys():
		l_unsorted = dHash[sz_hash]
		l_sorted = ivMergeSort(l_unsorted)
		dHash[sz_hash] = l_sorted

	return dHash

def hashGtfList(l_features,n_bin=16000):

	#
	# variables
	dHash = {}

	for lf in l_features:
		sz_hash = getFeatureHash(lf[_CHROM],lf[_START],n_bin=n_bin)
		if sz_hash not in dHash:
			dHash[sz_hash] = []

		dHash[sz_hash].append(copyFeature(lf))

	return dHash

def hashLookup(lf,dhash):

	# crete hash for feature's start and end
	shash = getFeatureHash(lf[_CHROM],lf[_START])
	ehash = getFeatureHash(lf[_CHROM],lf[_END])

	l_hits = []

	if shash not in dhash and ehash not in dhash:
		return set([])

	if shash == ehash:
		if shash in dhash:
			l_features = dhash[shash]
			n_features = len(l_features)
			i = 0
			while i < n_features and l_features[i][_START] < lf[_END]:
				if featureOverlap(l_features[i],lf):
					l_hits.append(l_features[i][_ATTR]['transcript_id'])

				i += 1

			sys.stderr.write("{} iterations\n".format(i))
	else:
		if shash in dhash:
			l_features = dhash[shash]
			n_features = len(l_features)
			i = 0
			while i < n_features and l_features[i][_START] < lf[_END]:
				if featureOverlap(l_features[i],lf):
					l_hits.append(l_features[i][_ATTR]['transcript_id'])

				i += 1

			sys.stderr.write("{} iterations\n".format(i))
		
		if ehash in dhash:
			l_features = dhash[ehash]
			n_features = len(l_features)
			i = 0
			while i < n_features and l_features[i][_START] < lf[_END]:
				if featureOverlap(l_features[i],lf):
					l_hits.append(l_features[i][_ATTR]['transcript_id'])

				i += 1

			sys.stderr.write("{} iterations\n".format(i))


	return set(l_hits)



#==============================================================================
# gtfToDict
#
# Parse a GTF file into a dict keyed by sz_key.
#==============================================================================
def gtfToDict(sz_file,sz_key="transcript_id"):

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
		szl = szl.strip()
		ll = szl.split("\t")

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
		lnew.append(ll[6]) # strand
		lnew.append(ll[2]) # type
		lnew.append(dattr.copy()) # attributes dict
		
		szId = dattr[sz_key]
		if szId not in d_out:
			d_out[szId] = lnew


	fin.close()

	return d_out
		

#==============================================================================
# gtfToList
#
# Parse a GTF file into a list with one entry per feature
#==============================================================================
def gtfToList(sz_file,sz_type="exon"):

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
		
#==============================================================================
# indexSortedGTFList
# Build a dict index of a sorted (by chrom and start position) list of 
# GTF features. Index holds start and end indices for each chromosome of 
# features in lgtf.
#==============================================================================
def indexSortedGTFList(lgtf):

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


#==============================================================================
# ivMergeSort
# Mergesort a bunch of parsed GTF features in a list
#==============================================================================
def ivMergeSort(lv):
	# if list is length 1 then it's sorted - sweet
	if len(lv) == 1:
		return lv

	# else we need to split the list in two and pass both halves to this function
	n_mid = len(lv)/2
	l_left = [lv[i] for i in range(n_mid)]
	l_right = [lv[i] for i in range(n_mid,len(lv))]

	# recurse
	l_left = ivMergeSort(l_left)
	l_right = ivMergeSort(l_right)

	# merge lists
	return ivMerge(l_left,l_right)

#==============================================================================
# ivMerge
# Merge segment of mergesort function
#==============================================================================
def ivMerge(ll,lr):
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



def copyFeature(lf):

	l_new = [lf[i] for i in range(len(lf)-1)]
	l_new.append(lf[_ATTR].copy())
	return l_new

def featureOverlap(lf1,lf2):

	if lf1[_START] < lf2[_END] and lf1[_END] > lf2[_START]:
		return True

	return False


def findOverlap(lfind,lcomp,lgtf,dindex,asSet=True,szId="transcript_id"):

	if lfind[_CHROM] not in dindex:
		raise Exception("invlaid chromosome request!")

	# jump to start of chromosome
	sz_chrom = lfind[_CHROM]

	n_start = dindex[sz_chrom][0]
	n_end = dindex[sz_chrom][1]
	
	# see if we can narrow the iteration window a little
	
	n_lastStart = n_start
	n_lastEnd = n_end
	b_done = False

	# if start is a hit we don't need to loop
	if featureOverlap(lcomp[n_start],lfind):
		b_done = True

	while not b_done:
		n_mid = (n_start+n_end)/2

		if lcomp[n_mid][_START] > lfind[_END]:
			# too far to the right, so go left
			# print "right",n_start,n_mid,n_end
			n_end = n_mid
		elif lcomp[n_mid][_START] < lfind[_END]:
			# start of feature is left of end of alignment, is it a hit?
			if lcomp[n_mid][_END] > lfind[_START]:
				# mid is a hit, is the start a hit too?
				if featureOverlap(lcomp[n_start],lfind):
					sys.stderr.write("warning: wtf encountered\n")
					n_start = n_lastStart
					
				#n_end = n_mid
				
				# mid is a hit so we're finished
				b_done = True
			else:
				# not a hit, go right
				n_start = n_mid

		if n_lastStart == n_start and n_lastEnd == n_end:
			# no more changes have been made so we're done - this is probably what happens
			# when there isn't an overlap
			b_done = True

		n_lastStart = n_start
		n_lastEnd = n_end

	#print i,n_end

	l_hits = []
	i = n_start
	ni = 0

	# loop 
	while i < dindex[sz_chrom][1] and lcomp[i][_START] < lfind[_END]:
		# check for overlap
		#if lcomp[i][_START] < lfind[_END] and lcomp[i][_END] > lfind[_START]:
			#l_hits.append(lcomp[i][_ATTR][szId])
			#print lgtf[i][_ATTR][szId]
		if featureOverlap(lcomp[i],lfind):
			#print lcomp[i]
			# l_comp points to a position in lgtf where the original elements live
			# so we have to scan that to see which of those elements are a hit
			ni = lcomp[i][-1]
			while lgtf[ni][_START] < lfind[_END]:
				#print lgtf[ni]
				if featureOverlap(lgtf[ni],lfind):
					l_hits.append(lgtf[ni][_ATTR][szId])

				ni += 1

		i += 1

	return set(l_hits)


#==============================================================================
# mergeList
# Merge overlapping adjacent features to produce a list of uniqe non overlapping
# features.
#==============================================================================
def mergeList(ll):

	l_out = []
	l_new = None
	n_index = 0

	for lf in ll:

		if l_new is None:
			l_new = copyFeature(lf)
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
				l_new = copyFeature(lf)
				l_new.append(n_index)

		n_index += 1

	# append final feature
	l_out.append(l_new)

	return l_out

def findOverlap_o2(lfind,lgtf,dindex,asSet=True,szId="transcript_id"):

	if lfind[_CHROM] not in dindex:
		raise Exception("invlaid chromosome request!")

	# jump to start of chromosome
	sz_chrom = lfind[_CHROM]

	n_start = dindex[sz_chrom][0]
	n_end = dindex[sz_chrom][1]
	
	# see if we can narrow the iteration window a little
	
	n_lastStart = n_start
	n_lastEnd = n_end
	b_done = False
	while not b_done:
		n_mid = (n_start+n_end)/2

		if lgtf[n_mid][_START] > lfind[_END]:
			# too far to the right, so go left
			# print "right",n_start,n_mid,n_end
			n_end = n_mid
		elif lgtf[n_mid][_START] < lfind[_END]:
			# start of feature is left of end of alignment, is it a hit?
			if lgtf[n_mid][_END] > lfind[_START]:
				# hit, go left
				n_end = n_mid
				b_done = True
			else:
				# not a hit, go right
				n_start = n_mid

		if n_lastStart == n_start and n_lastEnd == n_end:
			# no more changes have been made so we're done - this is probably what happens
			# when there isn't an overlap
			b_done = True

		n_lastStart = n_start
		n_lastEnd = n_end

	#print i,n_end

	l_hits = []
	i = n_start

	# loop 
	while i < dindex[sz_chrom][1] and lgtf[i][_START] <= lfind[_END]:
		# check for overlap
		if lgtf[i][_START] < lfind[_END] and lgtf[i][_END] > lfind[_START]:
			l_hits += lgtf[i][_ATTR][szId].split(",")
			#print lgtf[i][_ATTR][szId]

		i += 1

	return set(l_hits)


def findOverlap_oo(lfind,lgtf,dindex,asSet=True,szId="transcript_id"):

	if lfind[_CHROM] not in dindex:
		raise Exception("invlaid chromosome request!")

	# jump to start of chromosome
	sz_chrom = lfind[_CHROM]

	i = dindex[sz_chrom][0]
	n_end = dindex[sz_chrom][1]
	
	# see if we can narrow the iteration window a little
	
	b_done = False
	while not b_done:
		
		if n_end-i == 1:
			b_done = True

		else:
			
			n_mid = (n_end-i)/2 + i

			# overlap?
			if lgtf[n_mid][_START] < lfind[_END] and lgtf[n_mid][_END] > lfind[_START]:
				b_done = True
				# rewind until no overlap
				i = n_mid
				while lgtf[i][_START] < lfind[_END] and lgtf[i][_END] > lfind[_START]:
					i -= 1

			else:
				# no overlap - which way should we go?
				if lgtf[n_mid][_START] < lfind[_START]:
					# gtf is left so we go right
					#print "right",i,n_mid,n_end
					i = n_mid
				else:
					# gtf is right so go left
					#print "left",i,n_mid,n_end
					n_end = n_mid

	#print i,n_end

	l_hits = []

	# loop 
	while i < n_end and lgtf[i][_START] <= lfind[_END]:
		# check for overlap
		if lgtf[i][_START] < lfind[_END] and lgtf[i][_END] > lfind[_START]:
			l_hits += lgtf[i][_ATTR][szId].split(",")
			print lgtf[i][_ATTR][szId]

		i += 1

	return set(l_hits)

#n_ts = time.time()
#d_gtf = gtfToDict(sys.argv[1])
#n_te = time.time()
#print "Manual time: {:.4f} seconds".format(n_te-n_ts)

n_ts = time.time()
l_gtf = gtfToList(sys.argv[1])
n_te = time.time()
print "List time: {:.4f} seconds".format(n_te-n_ts)

n_ts = time.time()
l_gtfSorted = ivMergeSort(l_gtf)
n_te = time.time()
print "Sort time: {:.4f} seconds".format(n_te-n_ts)

#n_ts = time.time()
#d_hash = hashGtf(sys.argv[1])
#n_te = time.time()
#print "Hash time: {:.4f} seconds".format(n_te-n_ts)

n_ts = time.time()
d_hash = hashGtfList(l_gtfSorted)
n_te = time.time()
print "Hash time: {:.4f} seconds".format(n_te-n_ts)

l_lens = []
for l_bucket in d_hash:
	l_lens.append(len(l_bucket))

ar_lens = np.array(l_lens)

print ar_lens.mean(),ar_lens.min(),ar_lens.max()

#sys.exit()

# chr1:4,483,121-4,483,703
lfind = []
lfind.append(["chr1",6483278,6483378])
lfind.append(["chr1",106830754,106830856]) # 106,830,754-106,830,856
lfind.append(["chr1",4483238,4483338])
lfind.append(["chr1",4482500,4483338])

for lf in lfind:
	n_ts = time.time()
	s_hits = hashLookup(lf,d_hash)
	n_te = time.time()
	print "search completed in {:.8f}\n".format(n_te-n_ts),s_hits


