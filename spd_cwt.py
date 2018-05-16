#------------------------------------------------------------------------------
# spd_cwt.py
#
# Shawn Driscoll
# 20120328
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# CWT and support functions
#------------------------------------------------------------------------------

import numpy as np
import scipy.signal as ssig
from math import *

#------------------------------------------------------------------------------
# cwt_conv
#
# Creates a CWT spectrum via convolution
#------------------------------------------------------------------------------
def cwt_conv(x,srange=[1,10],steps=20,wavelet="morlet"):
	# set variables
	n_pnts = len(x)
	v_scales = np.zeros(steps)

	# scale increment - linear
	n_ds = (srange[1]-srange[0])/float(steps-1)

	if len(srange) < 2:
		return 1

	v_scales = np.arange(srange[0],srange[1],(srange[1]-srange[0])/float(steps-1))
	v_scales.resize(len(v_scales)+1)
	v_scales[-1] = v_scales[-2]+n_ds

	# create output matrix
	ar_out = np.array(np.zeros((steps,n_pnts)),dtype=complex)

	# create spectrum
	for i in range(steps):
		n_len = min(10*v_scales[i],n_pnts)
		if wavelet=="morlet":
			mo = ssig.morlet(v_scales[i],w=6)
		else:
			mo = ssig.morlet(v_scales[i],w=6)
		ar_out[i,:] = ssig.convolve(x,mo,mode='same')
	#done

	return ar_out

#end

