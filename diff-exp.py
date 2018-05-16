#!/usr/bin/env python
#
# diff-exp.py
#
# shawn driscoll
# 20120928
#
# differential expression on count data across two conditions with 
# or without biological replicates.
#

import sys,argparse
sys.path.append("/Users/pfafflab/coding/python")
import numpy as np
import bootstrap as bs
import de_utils
from math import log,sqrt,floor,ceil

#
# parse arguments
#

parser = argparse.ArgumentParser(description="Differential expression testing on count data across two conditions with or without replicates. Source files should have been produced by nexpr.")
parser.add_argument('condition_1_files',type=str,help="Comma separated list of files for condition 1")
parser.add_argument('condition_2_files',type=str,help="Comma separated list of files for condition 2")
parser.add_argument('-L',type=str,dest="labels",default="cond_1,cond_2",help="Conditions separated by comma (default: cond_1,cond_2)")
parser.add_argument('-k',type=float,dest='k_factor',default=0.5,help="Keep factor. The ratio of replicates from either condition that must have the minimum required hits. (default: 0.5)")
parser.add_argument('--sig-only',dest="sig_only",action="store_const",default=False,const=True,help="Return only significant genes (default: False)")

args = parser.parse_args()

#
# globals
#
MIN_COUNT = 10
NUM_BINS = 40

#
# main
#
def main(args):
	#
	# variables and initalizations
	#
	cond_1_files = args.condition_1_files.split(",")
	cond_2_files = args.condition_2_files.split(",")
	all_files = cond_1_files + cond_2_files

	num_cond_1 = len(cond_1_files)
	num_cond_2 = len(cond_2_files)
	num_samples = num_cond_1 + num_cond_2

	num_raw = 0
	num_masked = 0

	gene_info = []
	gene_info_masked = []
	gene_counts = None
	gene_counts_masked = None
	gene_counts_norm = None
	gene_expr = None
	gene_expr_masked = None

	row_mask = None

	i = 0
	j = 0
	k = 0

	condition_labels = args.labels.split(",")

	#
	# check file counts
	#
	if num_cond_1 == 0 or num_cond_2 == 0:
		sys.stderr.write("Error: missing input files!\n")
		return 1
	#
	# check condition label counts
	#
	if len(condition_labels) != 2:
		sys.stderr.write("Error: insufficient label count for conditions\n")
		return 1

	#
	# load files
	#
	sys.stderr.write("> loading files...\n")
	for i in range(num_samples):
		try:
			fp = open(all_files[i],"r")
		except IOError,e:
			sys.stderr.write("failed to open input file {:s}\n".format(all_files[i]))
			return 1

		# continue...
		k = 0
		for szl in fp:
			ll = szl.strip().split("\t")

			if i==0:
				#
				# first file! record info for this gene
				#
				gene_info.append([ll[0],ll[1],ll[2]])

				# build matrix
				if gene_counts is None:
					gene_counts = np.zeros(num_samples,dtype=float)
					gene_expr = np.zeros(num_samples,dtype=float)
				else:
					gene_counts = np.vstack((gene_counts,np.zeros(num_samples,dtype=float)))
					gene_expr = np.vstack((gene_expr,np.zeros(num_samples,dtype=float)))

			if len(np.shape(gene_counts)) == 1:
				gene_counts[i] = float(ll[3])
				gene_expr[i] = float(ll[4])
			else:
				gene_counts[k][i] = float(ll[3])
				gene_expr[k][i] = float(ll[4])

			k += 1

		# close file
		fp.close()

	#
	# files are loaded now we want to filter out low count genes by generating a mask
	#
	sys.stderr.write("> masking low count tags...\n")
	row_mask = condition_mask(gene_counts,num_cond_1,num_cond_2,args.k_factor)
	num_raw = np.shape(gene_counts)[0]
	num_masked = int(np.sum(row_mask))
	#
	# make the masked copies of the data
	#
	for i in range(num_raw):
		if row_mask[i] != 0:
			if gene_counts_masked is None:
				gene_counts_masked = gene_counts[i][:]
				gene_expr_masked = gene_expr[i][:]
			else:				
				gene_counts_masked = np.vstack((gene_counts_masked,gene_counts[i][:]))
				gene_expr_masked = np.vstack((gene_expr_masked,gene_expr[i][:]))

			gene_info_masked.append(gene_info[i])

	#
	# calculate normalization factors
	#
	norm_factors = np.zeros(num_samples)
	for i in range(num_samples):
		norm_factors[i] = float(np.percentile(gene_counts_masked[:,i],75))

	x = norm_factors.mean()
	norm_factors /= x

	#
	# normalize masked counts
	#
	for i in range(num_samples):
		gene_counts_masked[:,i] /= norm_factors[i]
	#
	# establish confidence interval bins
	#
	sys.stderr.write("> establishing confidence intervals...\n")

	# we'll divide up the data starting from the 1st percentile
	bin_bounds = np.percentile(gene_counts_masked,[1,99])
	bin_bounds[0] = transform_count(bin_bounds[0])
	bin_bounds[1] = transform_count(bin_bounds[1])

	bin_range = bin_bounds[1] - bin_bounds[0]
	bin_interval = bin_range/float(NUM_BINS)
	bin_matrix = np.zeros(NUM_BINS*2).reshape((NUM_BINS,2))

	for i in range(NUM_BINS):
		if i == 0:
			bin_matrix[i][0] = 0
		else:
			bin_matrix[i][0] = bin_bounds[0]+i*bin_interval

		bin_matrix[i][1] = bin_bounds[0]+(i+1)*bin_interval

	bin_matrix[NUM_BINS-1,1] = ceil(transform_count(gene_counts_masked.max()))

	diff_list = [[] for i in range(NUM_BINS)]	
	fold_list = [[] for i in range(NUM_BINS)]

	# generate distributions of absolute differences and absolute log2 fold changes
	# for all pairwise comparisons
	for i in range(num_masked):
		for j in range(num_samples-1):
			x1 = gene_counts_masked[i][j]
			b1 = int(hash_bin_index(transform_count(x1),bin_interval,bin_bounds[0],NUM_BINS-1))

			for k in range(j+1,num_samples):
				x2 = gene_counts_masked[i][k]
				b2 = int(hash_bin_index(transform_count(x2),bin_interval,bin_bounds[0],NUM_BINS-1))

				if x1 > 0 or x2 > 0:
					# difference
					diff = abs(x1-x2)
					diff_list[b1].append(diff)
					if b1 != b2:
						diff_list[b2].append(diff)

					# fold change
					fold = de_utils.log2_ratio(gene_counts_masked[i][j],gene_counts_masked[i][k])
					if np.isfinite(fold):
						fold = abs(fold)
						fold_list[b1].append(fold)
						if b1 != b2:
							fold_list[b2].append(fold)

	# lists are complete. now we can take the 95th percentile
	diff_lim = np.zeros(NUM_BINS)
	fold_lim = np.zeros(NUM_BINS)

	for i in range(NUM_BINS):
		if len(diff_list[i]) > 0:

			if num_samples > 2:
				temp = bs.bootstrap(diff_list[i],quar_stat_low)
				diff_lim[i] = temp.mean
			else:
				temp = bs.bootstrap(diff_list[i],quar_stat_high)
				diff_lim[i] = temp.mean


		if len(fold_list[i]) > 0:

			if num_samples > 2:
				temp = bs.bootstrap(fold_list[i],quar_stat_low)
				fold_lim[i] = temp.mean
			else:
				temp = bs.bootstrap(fold_list[i],quar_stat_high)
				fold_lim[i] = temp.mean

			# sys.stderr.write("{:f} {:f}\n".format(quar_stat_low(fold_list[i]),fold_lim[i]))


	#
	# do differential expression
	#
	sys.stderr.write("> running differential expression test...\n")

	cond_fold = np.zeros(num_masked)
	cond_sig = np.zeros(num_masked,dtype=int)
	cond_a_mean = np.zeros(num_masked)
	cond_b_mean = np.zeros(num_masked)

	for i in range(num_masked):
		# get mean count values
		cond_a_mean[i] = gene_counts_masked[i,0:num_cond_1].mean()
		cond_b_mean[i] = gene_counts_masked[i,num_cond_1:].mean()
		b1 = int(hash_bin_index(transform_count(cond_a_mean[i]),bin_interval,bin_bounds[0],NUM_BINS-1))
		b2 = int(hash_bin_index(transform_count(cond_b_mean[i]),bin_interval,bin_bounds[0],NUM_BINS-1))

		# do fold change
		cond_fold[i] = de_utils.log2_ratio(cond_a_mean[i],cond_b_mean[i])

		# get difference
		diff = abs(cond_b_mean[i]-cond_a_mean[i])

		# is this change tolerated?
		if np.isfinite(cond_fold[i]):
			#
			# fold change is finite
			#
			conf_bound_a = fold_lim[b1]
			conf_bound_b = fold_lim[b2]

			if abs(cond_fold[i]) > conf_bound_a and abs(cond_fold[i]) > conf_bound_b:
				cond_sig[i] = 2
			elif abs(cond_fold[i]) > conf_bound_a or abs(cond_fold[i]) > conf_bound_b:
				cond_sig[i] = 1
		else:
			#
			# fold change is infinite
			#
			conf_bound_a = diff_lim[b1]
			conf_bound_b = diff_lim[b2]

			if diff > conf_bound_a and diff > conf_bound_b:
				cond_sig[i] = 2
			elif diff > conf_bound_a or diff > conf_bound_b:
				cond_sig[i] = 1

	#
	# print out results!
	#
	sys.stdout.write("gene_id\ttranscript_id\tlocation\t{:s}\t{:s}\tlog2FoldChange\tsig_score\tsig\n".format(condition_labels[0],condition_labels[1]))
	for i in range(num_masked):
		if (args.sig_only and cond_sig[i] > 0) or not args.sig_only:
			lout = [
				"\t".join(gene_info_masked[i]),
				"{:0.4f}".format(cond_a_mean[i]),
				"{:0.4f}".format(cond_b_mean[i]),
				"{:0.4f}".format(cond_fold[i]),
				"{:d}".format(cond_sig[i])]

			if cond_sig[i] > 0:
				lout.append("yes")
			else:
				lout.append("no")

			sys.stdout.write("\t".join(lout) + "\n")


#
# condition_mask
# generates a set of transcript ids that will be kept for analysis. this 
# function checks each sample grouped by condition to determin if a 
# transcript has enough hits across samples per condition to be considered
# for analysis.
#
def condition_mask(data,num_cond_1,num_cond_2,k_factor):

	mask_out = []
	cond_1_values = []
	cond_2_values = []
	values = []
	num_rows = np.shape(data)[0]
	mask = np.zeros(num_rows)

	j = 0
	i = 0
	k = 0

	#
	# loop through data
	# 
	for i in range(num_rows):
		#
		# at each row we check counts across all samples
		cond_1_count = 0
		cond_2_count = 0
		
		# condition 1
		k = 0
		j = 0
		while k < num_cond_1:
			if data[i][k+j] >= MIN_COUNT:
				cond_1_count += 1
			k += 1
		# condition 2
		k = 0
		j = num_cond_1
		while k < num_cond_2:
			if data[i][k+j] >= MIN_COUNT:
				cond_2_count += 1
			k += 1

		if cond_1_count >= num_cond_1*k_factor or cond_2_count >= num_cond_2*k_factor:
			# passed, append key
			mask[i] = 1

	return mask

def transform_count(x):
	temp = x+1
	return log(temp)/log(10)

def hash_bin_index(x, iv, p25, max_bin):
	bin_index = floor( (x-p25)/iv )

	# check for bins beyond limit
	if bin_index > max_bin:
		return max_bin

	if bin_index < 0:
		return 0

	return bin_index

def hash_bin_index_list(x, iv, p25, max_bin):
	n = len(x)
	bin_index = np.zeros(n,dtype=int)
	i = 0

	for i in range(n):
		bin_index[i] = int(floor((x[i]-p25)/iv))
		if bin_index[i] > max_bin:
			bin_index[i] = max_bin
		elif bin_index[i] < 0:
			bin_index[i] = 0

	return bin_index



def diff_stat(x1,x2):
	return np.mean(x1)-np.mean(x2)

def quar_stat_low(x1):
	return np.percentile(np.array(x1),95)

def quar_stat_high(x1):
	return np.percentile(np.array(x1),97.5)


#
# entry point
#
if __name__ == "__main__":
	sys.exit(main(args))
