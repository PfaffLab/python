#
# de_utils.py
#
# shawn driscoll
# 20120928
#
# various functions useful for differential expression testing
#

import numpy as np
from math import log

#
# log2
# returns the log2 of value x
#
def log2(x):
	return log(x)/log(2)

#
# log2_ratio
# returns the log2 ratio of y/x
#
def log2_ratio(x,y):

	if x > 0 and y > 0:
		# x and y are both non-zero so we can do the normal ratio
		return log2(y/(x*1.0))
	else:
		# not both x and y are non-zero so check the conditions
		if y == 0 and x > 0:
			# y is zero
			return np.NINF
		elif x == 0 and y > 0:
			# x is zero
			return np.inf
		else:
			# both are zero
			return 0

	# return nan if something was missed
	return np.nan

