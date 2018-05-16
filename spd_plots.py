#==============================================================================
# spd_plots.py
#
# Shawn Driscoll
# 20120328
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Plotting functions to make plotting easier
#==============================================================================

import matplotlib.pyplot as plt

#==============================================================================
# lplot
# Quick line plot
#==============================================================================
def lplot(x):
	plt.plot(x)
	plt.show()


#==============================================================================
# lpplot
# plot first vector with the second on top of it plotted as dots
#==============================================================================
def lpplot(y1,y2,x=None):
	
	if x == None:
		# no x-axis values specified
		plt.plot(y1,'b-',y2,'ro')
		plt.show()
	else:
		plt.plot(x,y1,'b-',x,y2,'ro')
		plt.show()


#==============================================================================
# plotCwt
# Plot a cwt spectrum
#==============================================================================
def plotCwt(arCwt):
	plt.imshow(arCwt.real,aspect="auto")
	plt.show()


#==============================================================================
# plotCols
# Plot all columns of a matrix
#==============================================================================
def plotCols(ar):
	n = ar.shape[1]
	
	plt.clf()
	for i in range(n):
		plt.subplot(n,1,i+1)
		plt.plot(ar[:,i])

	plt.show()
	

#==============================================================================
# pairs
# Scatter plot all columns
#==============================================================================
def pairs(ar):
	n = ar.shape[1]
	npi = 1
	plt.clf()
	for i in range(n):
		for k in range(n):
			plt.subplot(n,n,npi)
			if k == i:
				plt.hist(ar[:,i])
			else:
				plt.plot(ar[:,k],ar[:,i],'bo')
			npi += 1

	plt.show()
	