#!/usr/bin/env python
#------------------------------------------------------------------------------
# tiffStackPca.py
#
# Shawn Driscoll
# 20120320
#
# Loads a tiff stack from file and performs PCA to find a set of dominant
# signals.
#------------------------------------------------------------------------------

import sys
import numpy as np
import Image as im
import ImageOps as iops
import mdp
from scipy.cluster.vq import vq, kmeans, whiten
import argparse
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# parse arguments
#------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Produce set of signals from a timeseries tiff stack using PCA and kmeans clustering")

parser.add_argument('stack', action="store", type=str, help="8 bit greyscale Tiff stack")
parser.add_argument('-k', action="store", dest="n_clusters", type=int, default=0, help="Number of clusters to use (default: auto)")
parser.add_argument('--standardize', action="store_const", dest="standardize", const=True, default=False, help="Standardize signal variances (default:False)")
parser.add_argument('-binx', action="store", dest="binx", type=int, default=1, help="bin this number of pixels in the x-dimension (default: 1)")
parser.add_argument('-biny', action="store", dest="biny", type=int, default=1, help="bin this number of pixels in the y-dimension (default: 1)")
parser.add_argument('-reshape', action="store", dest="sz_reshape", type=str, default="raw", help="bin image based on reshaping: ex: 128,128 (default: raw image is maintained)")

args = parser.parse_args()

#------------------------------------------------------------------------------
# variables
#------------------------------------------------------------------------------

b_calcBinFromReshape = False

if args.sz_reshape != "raw":
	ar_reshape = args.sz_reshape.split(",")
	b_calcBinFromReshape = True
#endif

# default max clusters - unless there is a user override there will never be more
# clusters used than this.
max_k = 20

ar_matrix = []
ar_trace = []

n = 0
n_frames = 0
n_traces = 0

sz_fbase = (args.stack).split(".")[0]
n_clusters = args.n_clusters

# create pca object
pca = mdp.nodes.PCANode(svd=True)

#------------------------------------------------------------------------------
# classes
#------------------------------------------------------------------------------

##
# stack iterator class
class ImageSequence:
	def __init__(self,im):
		self.im = im
	def __getitem__(self,ix):
		try:
			if ix:
				self.im.seek(ix)
			return self.im
		except EOFError:
			raise IndexError # end of sequence

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------

def binUnwrappedImage(ar,n_origx,n_origy,n_binx,n_biny):
	
	n_index=0
	n_binnedIndex=0

	# get dimensions
	m,n = np.shape(ar)
	# number of columns for binned image
	n_binnedColCount = n/n_binx/n_biny

	# make new matrix
	ar_new = np.zeros((m,n_binnedColCount))

	for j in range(n_origy):
		for i in range(n_origx):
			n_index = i+n_origx*j
			n_binnedIndex = int(np.floor(i/n_binx)+(n_origx/n_binx)*np.floor(j/n_biny))
			#print str(n_index) + "(" + str(j) + "," + str(i) + ") => " + str(n_binnedIndex)
			ar_new[:,n_binnedIndex] = ar_new[:,n_binnedIndex] + ar[:,n_index]
		#for
	#for

	ar_new = ar_new / (n_binx*n_biny)

	return ar_new

#def

def screeTest(ar_var):

	# create two new arrays
	ad1 = np.zeros(len(ar_var))
	ad2 = np.zeros(len(ar_var))

	# first derivative
	ad1[1:] = ar_var[1:] - ar_var[0:len(ar_var)-1]
	ad1[0] = ad1[1]
	# second derivative
	ad2[1:] = ad1[1:] - ad1[0:len(ad1)-1]
	ad2[0] = ad2[1]

	# find location of maximum
	n_maxLoc = ad2.argmax()

	return n_maxLoc+1
# end def

def varToPercentVar(ar_var):
	# create target array
	ar_percent = np.zeros(len(ar_var))
	# find sum of variances
	n_sum = sum(ar_var)
	# create new vector with each point of ar_var divided by the 
	# sum
	ar_percent[:] = ar_var[:]/n_sum
	# return the list
	return ar_percent
# end def

def estimateKmeansGroups(ar_data,maxg=15):
	# whiten data
	dwhite = whiten(ar_data)
	# get size of data array
	m,n = np.shape(ar_data)

	if maxg > m:
		maxg = m-1
	#if

	# initalize vector for sum of square results
	v_ss = np.zeros(maxg)

	# first grouping
	v_ss[0] = (m-1)*sum(np.var(dwhite,0))

	# iterate through clusterings to find variances at each partitioning
	for i in range(1,maxg):
		#np.random.seed((1000,3000))
		cents,vss = kmeans(dwhite,i+1)
		# find groupings
		code,dists = vq(dwhite,cents)
		# loop thorugh groupings
		for j in range(i+1):
			gmask = code != j
			# make mask
			arMask = np.ma.array(dwhite,mask=((gmask.repeat(n)).reshape(m,n)))
			n_cnt = arMask.count(0)[0]-1
			v_ss[i] = v_ss[i] + n_cnt*sum(np.var(arMask,0))
			#print str(j) + ", " + str(sum(gmask)) + ", " + str(v_ss[i])
		#endfor
		#print "\t".join(map(str,v_ss))	
	#for
	#print "\t".join(map(str,v_ss))
	return screeTest(v_ss)

# end def


#------------------------------------------------------------------------------
# main script
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# load stack into 2D array
#------------------------------------------------------------------------------

# open the stack
print "Loading stack from " + args.stack
imStack = im.open(args.stack)

# find dimension of the stack (height and x/y)
n_frames = 1
try:
	while 1:
		imStack.seek(imStack.tell()+1)
		n_frames = n_frames + 1
except EOFError:
	pass


# seek back to first frame, extract it and find the x, y dimensions of the stack.
imStack.seek(0)
ar_frame = np.array(imStack)
m_frame,n_frame = np.shape(ar_frame)
n_traces = len(ar_frame.flat)
print "> found " + str(n_frames) + " frames with " + str(n_traces) + " pixels per frame"

if b_calcBinFromReshape:
	args.binx = n_frame/(int(ar_reshape[1])*1.0)
	args.biny = m_frame/(int(ar_reshape[0])*1.0)
#endif

# confirm that binning divides the image evenly
if m_frame % args.biny != 0 or n_frame % args.binx != 0:
	print "*** Error ****\n> uneven binning. Image dimensions are " + str(m_frame) + " x " + str(n_frame) + ". Select a binning that divides the image evenly."
	sys.exit()
#endif

# create matrix
ar_matrix = np.zeros((n_frames,n_traces))
# rewind image
imStack.seek(0)

# loop through the frames and load values into 2D matrix. signals end up in 
# columns.

print "> building 2d matrix from 3d data..."
n = 0
for frame in ImageSequence(imStack):
	ar_frame = np.array(frame)
	# copy frame across columns at current now, n
	for i in range(n_traces):
		ar_matrix[n,i] = ar_frame.flat[i]
	# increment n
	n += 1

##
# bin the stack?

if args.binx > 1 or args.biny > 1:
	print "> binning pixels..."
	ar_pcaReady = binUnwrappedImage(ar_matrix,n_frame,m_frame,args.binx,args.biny)
else:
	# no binning, copy data to new matrix for pca prep
	ar_pcaReady = np.copy(ar_matrix)
#endif

##
# binned image dimensions?

m_frameBinned = int(m_frame/args.biny)
n_frameBinned = int(n_frame/args.binx)

#print np.shape(ar_matrix)
n_frames,n_tracesBinned = ar_pcaReady.shape
print "> traces reduced to " + str(n_tracesBinned)

#------------------------------------------------------------------------------
# PCA
#------------------------------------------------------------------------------

print "> prepping data..."
# zero mean the columns of the matrix 
# ** NOTE that this doesn't make a new copy of ar_matrix. in other words
# ar_matrix is actually modified in this process as ar_zm is only a "view"
# of ar_matrix.
ar_mean = np.mean(ar_pcaReady,0)  # this can also be used to create a mean image
ar_sd = np.std(ar_pcaReady,0)
m,n = np.shape(ar_pcaReady)

# subtract mean
#if args.b_noStand:
for i in range(n):
	ar_pcaReady[:,i] = ar_pcaReady[:,i] - ar_mean[i]
	# standardize?
	if args.standardize:
		ar_pcaReady[:, i] = ar_pcaReady[:, i] / ar_sd[i]

	#endfor
#else:
#	ar_pcaReady = whiten(ar_pcaReady)
#endif

# transpose
ar_pcaReady = ar_pcaReady.transpose()

print "> total variance = " + str(sum(np.var(ar_pcaReady,0)))

# start pca process
print "> pca..."
pca.train(ar_pcaReady)
# get projections
ar_scores = pca.execute(ar_pcaReady)
# square the eigenvalues
#pca.d = pow(pca.d,2)
# get eigenvectors
ar_eig = pca.v
#print np.shape(ar_eig)

##
# perform scree test to find number of significant components
n_sig = screeTest(pca.d)

##
# print the percent of variace explained by this number of significant components
ar_percent = varToPercentVar(pca.d)
print "> fraction of total variance for the first " + str(n_sig) + " principal components:"
print "> " + ", ".join(map(str,ar_percent[0:n_sig]))
print "> sum = " + str(sum(ar_percent[0:n_sig]))

#------------------------------------------------------------------------------
# kmeans on pca scores
#------------------------------------------------------------------------------

if args.n_clusters == 0:
	print "> estimating number of clusters for kmeans clustering..."
	n_clusters = estimateKmeansGroups(ar_scores[:,0:n_sig]) * 2
	
	if n_clusters > max_k:
		print "> warning: found " + str(n_clusters) + " clusters. defaulting to max of " + str(max_k)
		n_clusters = max_k
	else:
		print "> using " + str(n_clusters) + " clusters"

else:
	n_clusters = args.n_clusters

print "> running kmeans with " + str(n_clusters) + " clusters"
dwhite = whiten(ar_scores[:,0:n_sig])
cents,vss = kmeans(dwhite,n_clusters)
code,dists = vq(dwhite,cents)

print "> partitioning data based on clustering..."
# now partition the raw data based on the cluster groups
ar_means = np.zeros((n_frames,n_clusters))

for i in range(n_clusters):
	# break out group
	gmask = code != i
	arMask = np.ma.array(ar_pcaReady,mask=(gmask.repeat(n_frames)).reshape(n_tracesBinned,n_frames))
	t_mean = np.mean(arMask,0)
	ar_means[:,i] = t_mean[:]

#------------------------------------------------------------------------------
# saving results to disk
#------------------------------------------------------------------------------

print "> writing results to files..."

# reshape the code array and save it to disk
frameOut = code.reshape(np.floor(m_frame/args.biny),np.floor(n_frame/args.binx))
fout = open(sz_fbase + "_frame.txt",'w')
for i in range(int(np.floor(m_frame/args.biny))):
	fout.write("\t".join(map(str,frameOut[i,])) + "\n")
#for
fout.close()

# create mean image
imean = im.new('L',(m_frameBinned,n_frameBinned),"black")
imean.putdata(list(ar_mean))
imean = imean.convert('RGB')
#imean.save("mean_stack.tif")

# create images out of each cluster and plot as a grid of images
n_rows = n_clusters
n_rowsAdj = n_rows
n_cols = 1
# figure out layout such that there are more columns than rows
while n_cols < n_rowsAdj:
	n_cols = n_cols + 1
	n_rowsAdj = n_rows/n_cols

if n_rows % n_cols*n_rowsAdj > 0:
	#n_rowsAdj = n_rowsAdj+1
	n_cols = n_cols + 1

#fig = plt.figure()
for i in range(n_clusters):
	iout = im.new('L',(m_frameBinned,n_frameBinned),"black")
	g0 = code==i
	g0 = g0.astype(np.uint8)*255
	iout.putdata(list(g0))
	iout = iops.colorize(iout,(0,0,0),(255,0,0))
	#iout.save("cluster_0.tif")
	# blend mean of stack with this cluster's regions. flip image top to bottom
	ioutBlend = (im.blend(imean,iout,0.5)).transpose(im.FLIP_TOP_BOTTOM)
	#ioutBlend.save("cluster_" + str(i) + ".tif")
	#ax = fig.add_subplot(n_rowsAdj,n_cols,i+1,frame_on=False)
	ax = plt.subplot(n_rowsAdj,n_cols,i+1)
	plt.imshow(ioutBlend)
	plt.xlabel(str(i+1))
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.setp(ax.get_yticklabels(), visible=False)
plt.savefig(sz_fbase + "_clusters.pdf",dpi=300,format="pdf")
plt.clf()

# save traces out to text file
fout = open(sz_fbase + "_signals.txt",'w')
for i in range(n_frames):
	fout.write("\t".join(map(str,ar_means[i,])) + "\n")
#for
fout.close()

# create figure of the traces
#fig = plt.figure()
plt.subplots_adjust(hspace=0,wspace=0)
axl = plt.subplot(n_clusters,1,1,frame_on=True)
plt.plot(ar_means[:,0])
plt.ylabel("1")
plt.setp(axl.get_xticklabels(), visible=False)
plt.setp(axl.get_yticklabels(), visible=False)
for i in range(1,n_clusters):
	ax = plt.subplot(n_clusters,1,i+1,sharex=axl,frame_on=True)
	plt.plot(ar_means[:,i])
	plt.ylabel(str(i+1))
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.setp(ax.get_yticklabels(), visible=False)
#endfor
# make the last x-axis tick labels visible and small
plt.setp(ax.get_xticklabels(), visible=True)
plt.setp(ax.get_xticklabels(), fontsize=6)
plt.savefig(sz_fbase + "_signals.pdf",dpi=300,format="pdf")

# save data to file
#fout = open("pca_loadings.txt",'w')
#m,n = np.shape(ar_eig)
#for i in range(m):
#	fout.write("\t".join(map(str,ar_eig[i,0:n_sig])) + "\n")
#fout.close()

#fout = open("pca_scores.txt",'w')
#m,n = np.shape(ar_scores)
#for i in range(m):
#	fout.write("\t".join(map(str,ar_scores[i,0:n_sig])) + "\n")
#fout.close()

# find mean trace
#v_mean = np.mean(ar_zm,1)

#print "mean"
#for x in v_mean[:]:
#	print x

#print "mean"

##
# print the raw data out to STDOUT
#for n in range(n_frames):
	# print each row concatenated with tabs
#	print "\t".join(map(str,ar_matrix[n,]))

# ar = np.array(np.random.randint(256,size=(128,128)),dtype=np.uint8)

# put arra into an image
# img.putdata(list(ar.flat))
