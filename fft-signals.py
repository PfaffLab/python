#!/usr/bin/env python
#------------------------------------------------------------------------------
# fft-signals.py
#
# Shawn Driscoll
# 20120326
#
# Loads file passed at the command line and runs an FFT on each. The
# magnitude squared spectrum for each is plotted on a log scale.
#------------------------------------------------------------------------------

import sys,argparse
sys.path.append("/Users/pfafflab/coding/python")
import spd_signals as ss
import spd_utils as ssu
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
from math import *

#------------------------------------------------------------------------------
# parse arguments
#------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Produce magnitude squared FFT spectrum for each signal in input file. Signals in columns.")

parser.add_argument('signals', action="store", type=str, help="Tab delimited text file with signals in columns. Enter - to read from stdin.")
parser.add_argument('--sampling-rate', action="store", dest="n_samplingRate", type=float, default=1, help="Sampling rate in Hz. (default: 1 Hz)")
parser.add_argument('--skip-rows',action="store",dest="n_skipLines",type=int,default=0,help="Number of lines to skip from top of the file (default: 0)")
parser.add_argument('-range', action="store", dest="x_range", type=str, default="all", help="Specify range of frequencies to plot, eg 0.1,1 (default: all)")
parser.add_argument('--skip-cols', action="store", dest="n_firstCol", type=int, default=0, help="Number of columns to skip (eg skip a column of row labels) (default: 0)")
parser.add_argument("-o", action="store", dest="basename", type=str, default="", help="Basename for saved output file(s). Default: basename of input file or stdin")
parser.add_argument("--do-analysis", dest="doFFTAnalysis", action="store_const", const=True, default=False, help="Do FFT Analysis instead of producing power spectrum plots. (default: False)")

args = parser.parse_args()

#------------------------------------------------------------------------------
# variables
#------------------------------------------------------------------------------

ar_signals = []

if args.basename != "":
	sz_fbase = args.basename
else:
	if args.signals == "-":
		sz_fbase = "stdin"
	else:
		sz_fbase = (args.signals).split(".")[0]

ar_range = [0,args.n_samplingRate/2.0]
ar_pntRange = [0,0]

b_useRange = False

if args.x_range != "all":
	ar_range = map(float,args.x_range.split(","))
	b_useRange = True
#endif

#------------------------------------------------------------------------------
# main script
#------------------------------------------------------------------------------

# read data
sys.stderr.write("Loading signals from " + args.signals + "\n")
if args.signals == "-":
	ar_signals = np.loadtxt(sys.stdin,delimiter="\t",skiprows=int(args.n_skipLines))
else:
	ar_signals = np.loadtxt(args.signals,delimiter="\t",skiprows=int(args.n_skipLines))

n_len,n_signals = ar_signals.shape

if args.doFFTAnalysis:

	sys.stderr.write("> processing signals...\n")

	# make output array
	ar_out = np.array(np.zeros((n_signals-args.n_firstCol,3)))

	# loop through signals
	n_sig = 0
	for i in range(n_signals-args.n_firstCol):
		ar_res = ss.fftAnalysis(ar_signals[:,i+args.n_firstCol],samplingRate=args.n_samplingRate,v_range=ar_range)
		ar_out[n_sig,0] = i+args.n_firstCol
		ar_out[n_sig,1] = ar_res[0]
		ar_out[n_sig,2] = ar_res[1]
		n_sig += 1
	#done

	# save data to file
	ssu.writeTabDelim(sz_fbase + "_fftAnalysis.txt",ar_out)

else :

	n_lenPad = int(ss.nextPow2(n_len)/2 + 1)
	n_tpnts = int(ss.nextPow2(n_len)/2 + 1)
	n_step = (args.n_samplingRate/2.0) / (n_lenPad-1)
	v_freq = np.arange(n_tpnts,dtype=float)*n_step

	# if user specified a plot range in frequencies figure out the corresponding
	# point range to plot
	if b_useRange:
		ar_pntRange[0] = int(floor(ar_range[0]/n_step))
		ar_pntRange[1] = int(ceil(ar_range[1]/n_step))
	else:
		ar_range[0] = 0
		ar_range[1] = args.n_samplingRate/2.0
		ar_pntRange[0] = 0
		ar_pntRange[1] = n_lenPad
	# endif

	sys.stderr.write("> processing signals and creating plots...\n")

	# start plot
	plt.clf()
	plt.subplots_adjust(hspace=0,wspace=0)

	# loop through signals
	n_sig = 1
	for i in range(n_signals-args.n_firstCol):
		ar_ps,v_freq = ss.fftms(ar_signals[:,i+args.n_firstCol],samplingRate=args.n_samplingRate)

		# create red noise ps background
		rps = ss.redNoisePS(ar_signals[:,i+args.n_firstCol])
		# append plot
		ax = plt.subplot(n_signals-args.n_firstCol,1,n_sig)
		#plt.plot(ar_ps[ar_pntRange[0]:ar_pntRange[1],1],ar_ps[ar_pntRange[0]:ar_pntRange[1],0],'black',
		#	ar_ps[ar_pntRange[0]:ar_pntRange[1],1],rps[ar_pntRange[0]:ar_pntRange[1]],'r')
		plt.plot(v_freq[ar_pntRange[0]:ar_pntRange[1]],ar_ps[ar_pntRange[0]:ar_pntRange[1]],'black',
			v_freq[ar_pntRange[0]:ar_pntRange[1]],rps[ar_pntRange[0]:ar_pntRange[1]],'r')
		plt.xlim(ar_range[0],ar_range[1])
		plt.ylabel(str(n_sig))
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.setp(ax.get_yticklabels(),visible=False)

		n_sig += 1
	#done

	# activate the last x axis labels
	plt.setp(ax.get_xticklabels(),visible=True,fontsize=6)
	plt.xlabel("Frequency (cycles/unit-time)")
	plt.savefig(sz_fbase + "_fft.pdf",dpi=300,format="pdf")

