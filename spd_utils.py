#==============================================================================
# spd_utils.py
#
# Shawn Driscoll
# 20120327
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Convienence stuff
#==============================================================================

import numpy as np
from math import *
from scipy.stats import chi2

#==============================================================================
#==============================================================================
# FILE I/O
#==============================================================================
#==============================================================================


#==============================================================================
# writeTabDelim
# Exports an array as tab delim text. Can also handle 1D traces
#==============================================================================
def writeTabDelim(sz_file,x,header=None):
	
	n_dim = len(x.shape)
	n_rows = x.shape[0]

	# open output file
	fout = open(sz_file,'w')

	if header != None:
		# print header
		fout.write("\t".join(header) + "\n")

	if(n_dim == 1):
		# single column
		for i in range(n_rows):
			fout.write(str(x[i]) + "\n")
	else:
		for i in range(n_rows):
			# all columns in matrix
			fout.write("\t".join(map(str,x[i,:])) + "\n")

	fout.close()

#def



#==============================================================================
#==============================================================================
# Math Utilities
#==============================================================================
#==============================================================================

#==============================================================================
# findPeaksAboveSd
# use the windowed standard deviation as a detection method for finding 
# important peaks. works well for FFT spectrums.
#==============================================================================
def findPeaksAboveSd(x,nwin,k=4,prange=None):
	# variables
	s = x.std()
	nPnts = len(x)
	b_inPeak = False
	ar_peaks = []
	ar_positions = []
	nMax = 0
	nMaxPos = 0
	
	# get windowed standard deviation
	sdwin = windowSd(x,nwin)
	
	# if user passed range to prange use it, otherwise set it to the full range
	# of the data
	if prange == None:
		prange = [0,nPnts]
	
	# loop through the data
	for i in range(prange[1]-prange[0]):
		xi = i + prange[0]
		
		if sdwin[xi] > s*k:
			b_inPeak = True
			# is this position a peak or an inflected slope?
			if x[xi-1] < x[xi] and x[xi+1] < x[xi]:
				# it's a peak
				if x[xi] > nMax:
					nMax = x[xi]
					nMaxPos = xi
		else:
			if b_inPeak:
				b_inPeak = False
				if nMax > 0:
					ar_peaks.append(nMax)
					ar_positions.append(nMaxPos)
				nMax = 0
	
	return np.array(ar_peaks),np.array(ar_positions)

#==============================================================================
# windowSd
# Compute the standard deviation of an array in windows
#==============================================================================
def windowSd(x,nwin):
	# get length of data
	nPnts = len(x)
	
	# adjust window to the next odd value
	if nwin % 2. == 0:
		nwin += 1
	# half win for interations
	nwinHalf = np.floor(nwin/2)
	
	# variables for start and end of window at each iteration
	nwins = 0
	nwine = 0
	
	# copy x for output
	arOut = x.copy()
	
	for i in range(nPnts):
		nwins = i -nwinHalf
		nwine = nwins + nwin
		
		if nwins < 0:
			nwins = 0
		elif nwine > nPnts:
			nwine = nPnts
		
		arOut[i] = x[nwins:nwine].std()
	
	return arOut
	

#==============================================================================
# scale
# Like scale in R
#==============================================================================
def scale(ar,scale=True,center=True):
	m = ar.shape[0]
	
	# find mean and standard deviation of the columns
	rmean = ar.mean(0)
	rstd = ar.std(0)
	
	# make output array
	arOut = ar.copy()
	
	for i in range(m):
		if center:
			arOut[i] -= rmean
		if scale:
			arOut[i] /= rstd

	return arOut

#==============================================================================
# matMdist
# Returns a vector of mahalanobis distances of each vector in ars from the 
# center of mass of the vector set in arr.
# vectors are in rows (columns are dimensions)
#==============================================================================
def matMDist(ars,arr):
	
	# variables
	nSamples,nDim = ars.shape
	
	# create distance array
	mdist = np.arange(nSamples,dtype=float)
	
	# get mean vector of arr
	baseMean = arr.mean(0)	
	# subtract mean of arr from ars
	ars = np.subtract(ars,baseMean)
		
	# calculate covariance of arr
	arrzm = scale(arr,scale=False)
	acov = np.dot(arrzm.T,arrzm)
	acov /= arr.shape[0]
	# take determinant of covariance matrix (return this value)
	n_covDet = np.linalg.det(acov)
	
	# invert covariance matrix
	acov = np.linalg.inv(acov)
	
	# compute distance at each row of ars
	for i in range(nSamples):
		ap1 = np.dot(ars[i,:],acov)
		mdist[i] = np.dot(ap1,(ars[i,:]).T)
	
	return mdist,n_covDet

#==============================================================================
# mvoutlier
# Finds multivariate outliers using mahalanobis distance and, optionally, a subset of the
# original data (set robust=True)
#==============================================================================
def mvoutlier(ar,robust=False,alpha=0.05,inv=False,rinit=0.75):
	
	# dimension of data
	nSamples,nDim = ar.shape
	# subset dimension for robust analysis
	nh = int(nSamples*rinit)
	# threshold value
	nCutoff = chi2.isf(alpha/2,nDim)
	
	if robust:
		# perform robust analysis
		
		# subset data
		arh = np.random.permutation(ar)[0:nh]

		# get distance of all vectors from subset of vectors
		mdist,ndet = matMDist(ar,arh)
		ndetLast = ndet*2
		i = 0
		while ndet < ndetLast and i < 20:
			# re-subset the data based on the distances found. using
			# rinit*nSamples rows from those with the smallest distances
			# indices of distances
			indices = range(len(mdist))
			# sort indices on distance
			indices.sort(key = mdist.__getitem__)

			# create sorted copy of data
			arSorted = ar.copy()
			for i in range(nSamples):
				arSorted[i] = ar[indices[i]]
				
			# compute distances again based on new subset
			ndetLast = ndet
			mdist,ndet = matMDist(ar,arSorted[0:nh])
			
			i += 1
		
	else:	
		# get mahalanobis distance of each row vector from the mean row vector
		mdist,ndet = matMDist(ar,ar)
		
	# threshold
	
	if inv:
		# mask is true for the values that are NOT outliers
		return mdist <= nCutoff

	# mask is true for the outliers
	return mdist > nCutoff
		