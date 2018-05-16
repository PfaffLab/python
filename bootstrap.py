#
# bootstrap.py
#
# shawn driscoll
# 20120928
#
# bootstrap functions
#

import numpy as np
from math import sqrt

#
# objects
#
class BootstrapResult:
	mean = 0
	sem = 0
	conf_low = 0
	conf_high = 0
	data = []

	def __init__(self):
		self.mean = 0
		self.sem = 0
		self.conf_low = 0
		self.conf_high = 0
		self.it_results = []

	def stats(self):
		x = np.array(self.it_results)

		ptemp = np.percentile(x,[2.5,97.5])
		self.mean = np.mean(x)
		self.sem = sqrt(np.var(x))
		self.conf_low = ptemp[0]
		self.conf_high = ptemp[1]

#
# bootstrap
# bootstraps the mean of data over specified iterations
#
def bootstrap(data,stat_func,iterations=1000):

	# make sure data is a numpy array
	data_c = np.array(data)
	sample_size = len(data_c)
	bs_result = BootstrapResult()

	# iterate
	for i in range(iterations):
		indices = np.random.random_integers(0,sample_size-1,sample_size)
		temp = data_c[indices]
		bs_result.it_results.append(stat_func(temp))

	bs_result.stats()

	return bs_result

#
# two_sample_bootstrap_h0
#
# bootstraps the null hypothesis. for example the difference
# of the means of x1 and x2 is zero. stat_func would be a function
# that returns the difference between the mean of two arrays.
# the bootstrap distribution should be normal and centered at 
# zero. You can use the true difference of the means value
# within the bootstrap distribution to calculate a p-value.
#
def two_sample_bootstrap_h0(x1,x2,stat_func,iterations=10000):

	bs_r = BootstrapResult()
	i = 0
	length_1 = len(x1)
	length_2 = len(x2)
	length_full = length_1 + length_2

	ind_1 = []
	ind_2 = []

	x1c = np.array(x1)
	x2c = np.array(x2)

	# vector for result of iterations
	# bs_results = np.zeros(iterations)

	# create merged vector of all values
	all_values = np.append(x1,x2)

	# do it!
	for i in range(iterations):
		# 
		# create two random groups of values from the full set of values
		#
		ind_1 = random.sample(range(length_full),length_1)
		ind_2 = list(set(range(length_full)).difference(set(ind_1)))

		x1_p = all_values[ind_1]
		x2_p = all_values[ind_2]

		bs_r.it_results.append(stat_func(x1_p,x2_p))

	bs_r.stats()
	return bs_r


#
# two_sample_bootstrap
#
# for some function that operates on x1 and x2 and produces a 
# single number or statistic this function bootstraps said
# statistic for random samples of x1 and x2. if the resulting
# bootstrap distribution does not contain 0 in its 95% CI then
# it's likely the true result is valid.
#
def two_sample_bootstrap(x1,x2,stat_func,iterations=10000):

	bs_r = BootstrapResult()
	i = 0
	length_1 = len(x1)
	length_2 = len(x2)

	ind_1 = []
	ind_2 = []

	x1c = np.array(x1)
	x2c = np.array(x2)

	# vector for result of iterations
	# bs_results = np.zeros(iterations)

	# do it!
	for i in range(iterations):
		# 
		# sample with replacement from each group
		#
		x1_s = x1c[np.random.random_integers(0,length_1-1,length_1)]
		x2_s = x2c[np.random.random_integers(0,length_2-1,length_2)]

		bs_r.it_results.append(stat_func(x1_s,x2_s))

	bs_r.stats()
	return bs_r

