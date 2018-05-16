#==============================================================================
# spd_signals.py
#
# Shawn Driscoll
# 20120327
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Library of signal processing functions
#==============================================================================

import sys
sys.path.append("/Users/pfafflab/coding/python")
import spd_utils as su
import numpy as np
from numpy.fft import *
from random import gauss
import scipy.signal as ssig
from math import *
import matplotlib.pyplot as plt

#==============================================================================
# nextPow2
# Returns the number that is the next power of two past the value, n
#==============================================================================
def nextPow2(n):
	n_next = ceil(log(n,2))
	return 2 ** n_next
#def

#==============================================================================
# getPaddedSignal
# Returns a copy of the passed signal padded with zeros out to the next
# power of 2 points in length.
#==============================================================================
def getPaddedSignal(v_signal):

	n_pad = nextPow2(len(v_signal))
	ar_pad = np.array(np.zeros(n_pad),dtype=float)

	ar_pad[0:len(v_signal)] = v_signal[:]

	return ar_pad

# end def


#==============================================================================
# fftms
# Returns the magnitude squared power spectrum for the passed signal. signal
# is padded with zeros to the next largest power of two.
#==============================================================================
def fftms(v_signal,samplingRate=1,plot=False):

	# subtract mean
	n_mean = v_signal.mean()
	v_signal -= n_mean
	v_padded = getPaddedSignal(v_signal)
	n_tpoints = len(v_padded)/2 + 1
	n_step = (samplingRate/2.0)/(n_tpoints-1)

	# do fft
	sp = fft(v_padded)

	# create magnitude squared spectrum
	v_magSqr = sp.real ** 2 + sp.imag ** 2

	# return N/2 +1 points
	v_result = v_magSqr[0:n_tpoints]
	#v_result[0:] = v_magSqr[0:n_tpoints]
	v_result[0] = 0

	# normalize
	v_result /= len(v_signal)

	# create scaling vector
	v_freq = np.arange(n_tpoints,dtype=float) * n_step

	if(plot):
		plt.plot(v_freq,v_result)
		plt.show()
	#endif

	return v_result,v_freq

# end def


#==============================================================================
# autocovariance
# Returns the autocovariance of a signal at a specified lag
#==============================================================================
def autocovariance(v_signal,n_lag=1):
	# 1/N S (x(t) -- xbar)(x(t + m) -- xbar)
	# N = numpnts(w)
	# m = lag
	# x(t) = w
	# xbar = mean of w

	n_len = len(v_signal)
	n_mean = np.mean(v_signal)

	# create two new vectors
	v1 = np.zeros(n_len-n_lag)
	v2 = v1.copy()

	# copy data and lagged data
	v1[:] = v_signal[0:len(v1)]-n_mean
	v2[:] = v_signal[n_lag:]-n_mean

	v3 = v1 * v2

	return sum(v3)/n_len
#done

#==============================================================================
# acc
# Returns the lag N autocorrelation coefficient for the signal in v_signal
#==============================================================================
def acc(v_signal,n_lag=1):
	n_autoCov = autocovariance(v_signal,n_lag=n_lag)
	n_var = np.var(v_signal)
	return n_autoCov/n_var
#done

#==============================================================================
# redNoisePS
# Returns np.array containing the mean red-noise power spectrum
# for the signal, v_signal with lag-N correlation coefficient calculated
# for n_lag.
#==============================================================================
def redNoisePS(v_signal,n_lag=1):
	# get autocorrelation coefficient
	n_acc = acc(v_signal,n_lag=n_lag)
	n_pnts = len(v_signal)
	n_pad = int(nextPow2(n_pnts))
	n_tpnts = int(n_pad/2 + 1)
	n_var = np.var(v_signal)
	n_step = 0.5/(n_tpnts-1)

	# power spectrum will have nextPow2/2 +1 points
	v_ps = np.zeros(n_tpnts)
	v_freq = np.arange(n_tpnts,dtype=float) * n_step

	#done

	# wout[] = (1-nacc^2) / (1 + nacc^2 - 2*nacc*cos(2*pi*wout[p]))
	for i in range(n_tpnts):
		v_ps[i] = (1-n_acc**2) / (1 + n_acc**2 - 2*n_acc*cos(2*pi*v_freq[i]))
	#done

	# scale by data's variance
	v_ps *= n_var
	return v_ps
#done


#==============================================================================
# redNoiseTS
# Returns a red-noise timeseries based on the supplied signal with an acc
# at the supplied lag value
#==============================================================================
def redNoiseTS(v_signal,n_lag=1):
	# get autocorrelation coefficient
	n_acc = acc(v_signal,n_lag=n_lag)
	# solve v(ar1) = v(white)/(1-n_acc^2) where v(ar1) = v(v_signal)
	n_var = np.var(v_signal)
	n_whiteVar = n_var * (1-n_acc**2)

	# create output
	v_out = v_signal.copy()
	v_out[:] = 0

	for i in range(1,len(v_signal)):
		v_out[i] = n_acc * v_out[i-1] + gauss(0,sqrt(n_whiteVar))
	#done

	n_mean = np.mean(v_out)
	v_out[:] -= n_mean

	return v_out
#done

#==============================================================================
# fftAnalysis
# Creates magnitude squared FFT and matching red noise power spectrum 
# background as well as a background normalized copy of the power sepctrum.
#==============================================================================
def fftAnalysis(v_signal,samplingRate=1,n_lag=1,plot=False,v_range=[]):
	# set variables

	b_useRange = False
	if(len(v_range) > 0):
		b_useRange = True
	#endif

	# do fft
	ps,freq = fftms(v_signal,samplingRate=samplingRate)
	# create red noise background
	rps = redNoisePS(v_signal,n_lag=n_lag)

	# create red noise normalized copy of the power spectrum
	psMod = ps/rps

	# get maxima of psMod
	psModMaxima = maxima(psMod)

	# find position of highest frequency peak
	if b_useRange:
		# get frequency step interval
		n_step = freq[1]-freq[0]
		# convert frequency range into indexes
		v_rangei = [int(floor(v_range[0]/n_step)),int(ceil(v_range[1]/n_step))]
		n_peakIndex = index_max(psModMaxima[v_rangei[0]:v_rangei[1]]) + v_rangei[0]
	else:
		n_peakIndex = index_max(psModMaxima)
	#endif

	n_peakFreq = freq[n_peakIndex]
	n_peakScore = psMod[n_peakIndex]

	if plot:
		# create timing vector for time series
		v_samples = np.arange(len(v_signal),dtype=float)/samplingRate
		#endfor
		plt.subplot(311)
		plt.plot(v_samples,v_signal)
		plt.ylabel("Response")
		ax = plt.subplot(312)
		plt.plot(freq,ps,'black',freq,rps,'r')
		plt.xscale('log')
		plt.ylabel("Power")
		plt.subplot(313,sharex=ax)
		plt.plot(freq,psMod)
		plt.ylabel("Normalized Power")
		plt.show()
	#endif

	return n_peakFreq,n_peakScore

#done


#==============================================================================
# diff
# Take first time derivative of signal, x
#==============================================================================
def diff(x):

	xDiff = x.copy()
	n_pnts = len(x)

	for i in range(1,n_pnts):
		xDiff[i] = x[i]-x[i-1]
	#done

	xDiff[0] = xDiff[1]

	return(xDiff)
#def


#==============================================================================
# maxima
# Returns a copy of signal, x, with only maxima remaining. all other poitns
# set to nan
#==============================================================================
def maxima(x):
	# set variables
	n_pnts = len(x)

	# diff
	xD = diff(x)

	# create maxima vector and set all values to nan
	xMax = xD.copy()
	xMax[:] = float('nan')

	# find zero crossings in the derivative
	for i in range(1,n_pnts):

		if xD[i] == 0:
			if xD[i-1] > 0 and xD[i] < 0:
				xMax[i] = x[i]
			#endif
		else:
			if xD[i-1] > 0 and xD[i] < 0:
				xMax[i-1] = x[i-1]
			#endif
		#endif

	#done
	return(xMax)
#def

#==============================================================================
# index_max
# Returns the index of the point with the greatest value in x. Nan aware.
#==============================================================================
def index_max(x):
	# set variables
	n_pnts = len(x)
	n_maxi = -1
	n_max = 0

	for i in range(n_pnts):
		if not isnan(x[i]):
			if x[i] > n_max:
				n_maxi = i
				n_max = x[i]
		#endif
	#done

	return(n_maxi)
#def


#==============================================================================
# cwt
# Creates a CWT spectrum from signal. Specify scales as widths of the 
# wavelet (or periods of signal). s0 is the shorter period and sN is the
# last period (eg. 1,10 for periods 1s to 10s or 1 unit to 10 units).
# CWT resolution is controlled by numScales. Set power=True to return the 
# magnitude squared spectrum. qthresh can be used to apply quantile thresholding.
#==============================================================================
def cwt(signal,s0,sN,numScales,wavelet="morlet",power=False,scaling="linear",qthresh=False,samplingRate=1,plot=False):
	# variables
	npnts = len(signal)
	npad = nextPow2(npnts)
	# wavelet param
	param = 6
	# scale step
	dj = 0
	# sampling rate
	dt = 1./samplingRate

	# get power of 2 points padded signal and remove mean
	psignal = getPaddedSignal(signal) - signal.mean()

	# set param based on wavelet
	if wavelet=="morlet":
		param = 6
	elif wavelet == "paul":
		param = 4
	elif wavelet == "dog":
		param = 2
	else:
		wavelet = "morlet"
		param = 6
	#endif

	# setup scaling vector
	if scaling == "pow2":
		dj = log(sN/s0) / (numScales*log(2))
		scales = np.arange(numScales+1,dtype=float)
		scales = s0*2**(scales*dj)
	else:
		# linear scaling
		dj = (sN-s0)/float(numScales)
		scales = s0 + np.arange(numScales+1,dtype=float)*dj
	#endif

	#print dt
	#scales *= samplingRate

	# setup k vector

	# establish first 1/2 of k which is a linear ramp
	k = np.arange(npad/2+1,dtype=float)
	k = k*((2*pi)/(npad*dt))
	kn = k.copy()
	kn -= kn.max()
	# final k is concatenation of k and kn for a total of npad points
	k = np.concatenate((k,kn[1:npad/2]))

	#return k

	# do fft of the padded signal - python returns the symmetrical specturm npad points
	# long
	f = fft(psignal)

	# create empty period and wavelet arrays
	period = scales.copy()
	cwtOut = np.array(np.zeros((numScales+1,npad)),dtype=complex)

	# loop through scales and compute the transform
	for i in range(numScales+1):
		daughter,ffactor = wave_bases(wavelet,k,scales[i],param)
		cwtOut[i,:] = ifft(f * daughter)
	#done

	# adjust periods
	period = ffactor * scales

	# remove padding from cwt matrix
	cwtOut = cwtOut[:,0:npnts]

	if power:
		cwtOut = cwtOut.real**2 + cwtOut.imag**2
	#endif

	if plot:
		plt.clf()
		plt.subplot(311)
		plt.imshow(np.real(cwtOut),aspect="auto")
		plt.subplot(312)
		plt.plot(period,np.real(cwtOut).mean(1))
		plt.subplot(313)
		plt.plot(signal)
		plt.show()
	#endif

	return cwtOut,scales,period

#end

#==============================================================================
# wave_bases
# Support for CWT.  Creates "daughter" that is multiplied by the FFT spectrum
# at each scale of the CWT.
#==============================================================================
def wave_bases(wavelet, wk, scale, param):
	####
	# variables
	i=0
	n_pnts = len(wk)
	n_ffactor = 0

	# create vectors
	ar_expnt = wk.copy()



	if wavelet == "morlet":
		# fourier factor
		n_ffactor = (4*pi)/(param+sqrt(2+param**2))
		# expnt array
		ar_expnt = -(scale*wk - param)**2/2 * (wk > 0)
		# normalization
		nnorm = sqrt(scale*wk[1]) * (pi**(-0.25)) * sqrt(n_pnts)
		# create daughter
		ar_daughter = wk.copy()
		for i in range(n_pnts):
			ar_daughter[i] = nnorm * exp(ar_expnt[i])
		#done
		ar_daughter *= (wk > 0)
	elif wavelet == "paul":
		n_ffactor = 4*pi/(2*param+1)
		# expnt vector
		ar_expnt = -(scale*wk) * (wk > 0)

		nprod = 2
		i = 3
		while i < 2*param:
			nprod *= i
			i += 1

		# normalization
		nnorm = sqrt(scale*wk[1])*((2**param)/sqrt(param*nprod))*sqrt(n_pnts)

		# daughter
		ar_daughter = wk.copy()
		for i in range(n_pnts):
			ar_daughter[i] = nnorm*((scale*wk[i])**param)*exp(ar_expnt[i])
		#done
		ar_daughter *= (wk > 0)

	elif wavelet == "dog":
		# fourier factor
		n_ffactor = 2*pi*sqrt(2./(2*param+1))

		# expnt vector
		ar_expnt = -(scale*wk)**2 / 2.0
		# normalization
		nnorm = sqrt(scale*wk[1]/gamma(param+0.5))*sqrt(n_pnts)
		
		# create daughter
		ar_daughter = np.array(np.zeros(wk.shape),dtype=complex)
		for i in range(n_pnts):
			ar_daughter[i] = -nnorm*(1j**param)*((scale*wk[i])**param)*exp(ar_expnt[i])
		#done

	return ar_daughter,n_ffactor

#end

#==============================================================================
# morletToFourier
# convert from morlet scaling and fourier scaling
#==============================================================================
def morletToFourier(n,param=6):
	return (4*pi*n)/(param + sqrt(2 + param**2))
#def

#==============================================================================
# fourierToMorlet
# convert from fourier scaling to morlet scaling
#==============================================================================
def fourierToMorlet(n,param=6):
	return (n*(param + sqrt(2 + param**2)))/(4*pi)
#def

#==============================================================================
# cwtProperties
#
# find properties of the CWT spectrum (mean spectrum, ridgeline, etc)
#==============================================================================
def cwtProperties(arcwt):

	# is the cwt complex?
	if arcwt.dtype == 'complex':
		arcwt = arcwt.real**2 + arcwt.imag**2
	#endif

	# get dimensions
	scales,npnts = arcwt.shape

	# create mean spectrum
	cmean = arcwt.mean(1)
	# find peaks in the mean spectrum
	cmaxima = maxima(cmean)
	# peak scale
	n_peakScale = index_max(cmaxima)

	# peak ridge line
	ridge = np.arange(npnts,dtype=float)
	ridge[:] = float('nan')
	ra = ridge.copy()

	###
	# create time position peaks vector (ridge line)
	for i in range(npnts):
		# find peaks at current time slice
		cmaxima = maxima(arcwt[:,i])

		# find peak closest to the mean peak position
		n_minLoc = -1
		n_minDiff = scales
		n_diff = 0
		for n in range(scales):
			if not isnan(cmaxima[n]):
				n_diff = abs(n-n_peakScale)
				if n_diff < n_minDiff:
					n_minDiff = n_diff
					n_minLoc = n
				#endif
			#endif
			nPeak = n_minLoc
		#endfor

		if nPeak > 0:
			ridge[i] = nPeak
	#endfor

	ridgePower = ridge.copy()
	for i in range(npnts):
		if not isnan(ridge[i]):
			ridgePower[i] = arcwt[ridge[i],i]
		else:
			ridgePower[i] = float('nan')
		#endif
	#endfor

	return cmean,ridge,ridgePower

#def

#==============================================================================
# hpfilter
#
# IIR high pass filter applied forward and backwards to avoid phase drift
#==============================================================================
def hpfilter(signal,ff,f,offset=0.05,ftype='ellip'):

	# normalize frequency by nyquist frequency
	wp = ff*2./f
	ws = (ff-offset)*2./f

	# design filter
	b,a = ssig.iirdesign(wp,ws,gpass=1,gstop=10,ftype=ftype)

	# apply filter to signal
	sf = ssig.filtfilt(b,a,signal)

	# return filtered signal
	return sf
#end

#==============================================================================
# lpfilter
#
# IIR high pass filter applied forward and backwards to avoid phase drift.
# Parameters:
# ff  - cutoff frequency
# f   - sampling rate in Hz
#==============================================================================
def lpfilter(signal,ff,f,offset=0.05,ftype='ellip'):

	# normalize frequency by nyquist frequency
	wp = ff*2./f
	ws = (ff+offset)*2./f

	# design filter
	b,a = ssig.iirdesign(wp,ws,gpass=1,gstop=10,ftype=ftype)

	# apply filter to signal
	sf = ssig.filtfilt(b,a,signal)

	# return filtered signal
	return sf
#end

#==============================================================================
# bpfilter
#
# IIR band pass filter applied forward and backwards to avoid phase drift
# Parameters:
# f0   - start of passband (Hz)
# fn   - end of passband (Hz)
# f    - sampling rate of signal (Hz)
#==============================================================================
def bpfilter(signal,f0,fn,f,offset=0.05,ftype='ellip'):

	# normalize frequency by nyquist frequency
	wp = [f0*2./f,fn*2./f]
	ws = [(f0-offset)*2./f, (fn+offset)*2./f]

	# design filter
	b,a = ssig.iirdesign(wp,ws,gpass=1,gstop=10,ftype=ftype)

	# apply filter to signal
	sf = ssig.filtfilt(b,a,signal)

	# return filtered signal
	return sf
#end


#==============================================================================
# fcwts
# 
# This function looks for "signal" in timeseries data using a series of
# tests.
#==============================================================================
def fcwts(signal,s0,prange=[2,6]):

	dt = 1./s0
	fdt = 0
	nPeaks = 0
	zPeak = 0
	nfPeriod = 0
	
	b_fftPeak = False
	b_cwtPeriod = False
	b_cwtMag = False

	# make a copy of the input and detrend it with a hpfilter
	x = signal.copy()
	nMean = x.mean()
	x -= nMean
	x = hpfilter(x,0.1,s0)
	
	# adjust range for CWT and FFT
	nRangeFactor = (prange[1]-prange[0])*0.25
	crange = [prange[0]-nRangeFactor,prange[1]+nRangeFactor]
	
	# expected stats
	nExpectedMean = sum(prange)/2.
	nExpectedSdev = (prange[1]-nExpectedMean)/1.96
	
	# do FFT of signal to find out if there are any interesting peaks within 
	# the expected range
	ps,freq = fftms(x,samplingRate=s0)
	fdt = freq[1]-freq[0]
	lPeaks,lPos = su.findPeaksAboveSd(ps,np.floor(0.05/fdt))
	# figure out points corresponding to frequency boundaries of interest
	frange = [1./crange[1],1./crange[0]]
	fprange = [0,0]
	i = 0
	while freq[i] < frange[0]:
		i += 1
	fprange[0] = i
	while freq[i] < frange[1]:
		i += 1
	fprange[1] = i
	
	# examin detected peaks
	if len(lPeaks) > 0:
		nPeaks = len(lPeaks)
		if nPeaks == 1:
			# just one peak - is it good?
			nfPeriod = 1./freq[lPos[0]]
			zPeak = (nfPeriod-nExpectedMean)/nExpectedSdev
			if abs(zPeak) < 2:
				b_fftPeak = True
				# adjust expected range to center on this value
				
				
				

#==============================================================================
#==============================================================================
# BURST ANALYSIS FUNCTIONS
#==============================================================================
#==============================================================================

#==============================================================================
# ephysMakeBpm
# Make a burst position matrix for a signal. 
#==============================================================================
def ephysMakeBpm(signal,fs=1,prange=[2,6],autop=False,mthresh=False):
	
	# setup variables
	dt = 1/fs
	nScales = 50
	nPnts = len(signal)
	
	# bandpass filter the input signal
	fsignal = bpfilter(signal,0.1,1,fs)
	
	# do cwt of input singal
	arCwt,scales,periods = cwt(signal,prange[0],prange[1],nScales,samplingRate=fs)
	
	# get properties of the CWT
	cmean,cridge,cridgePower = cwtProperties(arCwt)
	
	# trace out ridge in the phase data of the CWT
	cphase = cridge.copy()
	for i in range(nPnts):
		if not isnan(cridge[i]):
			cphase[i] = atan2(np.imag(arCwt[cridge[i],i]),np.real(arCwt[cridge[i],i]))
		else:
			cphase[i] = float('nan')
		#endif
	#endfor
	
	# the phase information will be a sawtooth kinda wave but we have to figure
	# out which direction it's sloping. basically we want the location of the 
	# peaks in the sawtooth waveform
	cphasePeaks = cphase.copy()
	cphasePeaks[:] = -1
	# find slopes throughout the phase data
	cdiff = diff(cphase)
	
	if np.median(cdiff) > 0:
		# sawtooth with upward ramps (positive median slope)
		for i in range(nPnts-1):
			if cphase[i] > 0 and cphase[i+1] < 0:
				cphasePeaks[i] = i
	else:
		# sawtooth with downward ramps
		for i in range(1,nPnts):
			if cphase[i] > 0 and cphase[i-1] < 0:
				cphasePeaks[i] = i

	# find peak count
	cpCount = sum(cphasePeaks > 0)
	# create array to hold indices of the peaks
	cpi = np.array(np.zeros(cpCount),dtype=int)
	# put indicies in new array
	n = 0
	for i in range(nPnts):
		if cphasePeaks[i] > 0:
			cpi[n] = i
			n += 1
	
	# find peaks in the filtered data between these peaks
	cpeaks = np.arange(cpCount-1,dtype=int)
	for i in range(cpCount-1):
		# max location between two phase peaks
		npos = fsignal[cpi[i]:cpi[i+1]].argmax()
		cpeaks[i] = cpi[i]+npos
	
	# now with the location of the peaks we can find the minima between each peak
	cminima = np.arange(cpCount-2,dtype=int)
	for i in range(cpCount-2):
		# min location between two peaks
		npos = fsignal[cpeaks[i]:cpeaks[i+1]].argmin()
		cminima[i] = npos+cpeaks[i]
	
	# number of bursts will be one less than the number of minima
	bpm = np.array(np.zeros((len(cminima)-1,8)),dtype=float)
	
	# load up the matrix
	for i in range(len(cminima)-1):
		# start
		bpm[i,0] = cminima[i]
		# end
		bpm[i,1] = cminima[i+1]
		# peak
		bpm[i,2] = cpeaks[i+1]
		# amplitude
		bpm[i,3] = fsignal[cpeaks[i+1]] - fsignal[cminima[i]]
		# duration
		bpm[i,4] = bpm[i,1]-bpm[i,0]
		# find onset @ 25% of front amplitude
		dy = bpm[i,3]*0.25
		k = cminima[i]
		while fsignal[k] < fsignal[cminima[i]]+dy:
			k+=1
		bpm[i,5] = k
		# find offset @ 25% of back amplitude
		dy = (fsignal[cpeaks[i+1]]-fsignal[cminima[i+1]])*0.25
		k = cpeaks[i+1]
		while fsignal[k] > fsignal[cminima[i+1]]+dy:
			k += 1
		bpm[i,6] = k
		# burst width (offset position - onset position)
		bpm[i,7] = bpm[i,6]-bpm[i,5]
	
	# multidimensional thresholding of the amplitude/duration data to remove
	# extreme outliers. This protects the burst selection from funky data that
	# was interpreted as a cycle by the CWT but is in fact some type of junk
	if mthresh:
		# find outliers
		omask = su.mvoutlier(bpm[:,3:5],robust=True)
			
		if sum(omask) > 0:
			print "> dropping " + str(sum(omask)) + " rows from position matrix"
			# subset bpm by dropping outlier rows
			nKeep = sum(omask==False)
			bpmi = np.array(np.zeros((nKeep,bpm.shape[1])))
			n = 0
			for i in range(len(omask)):
				if not omask[i]:
					bpmi[n] = bpm[i]
					n += 1
			bpm = bpmi.copy()
	
	return bpm,["start","end","peak","amplitude","cycle_duration","onset","offset","burst_duration"]
	
#def

		
	
	
	