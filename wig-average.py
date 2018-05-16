#!/usr/bin/env python
#
# wig-average.py
#
# Shawn Driscoll
# 20120912
# 
# 

import sys,os
from math import ceil,floor,log
from struct import pack,unpack
import numpy as np
# import tempfile

NOISE_SHIFT = 1

def main(args):

	if len(args) < 3:
		sys.stderr.write("error: too few arguments!\n")
		return 1

	read_length = int(args[0])
	wig1 = args[1]
	wig2 = args[2]
	fp1 = None
	fp2 = None
	wig1_index = {}
	wig2_index = {}
	wig1_coverage = 0
	wig2_coverage = 0

	#
	# index the wigs
	#

	sys.stderr.write("indexing {:s}...\n".format(wig1))
	wig1_index,wig1_coverage = index_wig(wig1)
	wig1_coverage = wig1_coverage/(read_length*1.0)/1000000.0
	sys.stderr.write("> read coverage {:g} (millions)\n".format(wig1_coverage))

	sys.stderr.write("indexing {:s}...\n".format(wig2))
	wig2_index,wig2_coverage = index_wig(wig2)
	wig2_coverage = wig2_coverage/(read_length*1.0)/1000000.0
	sys.stderr.write("> read coverage {:g} (millions)\n".format(wig2_coverage))

	if wig1_index == None or wig2_index == None:
		sys.stderr.write("failed to index wig files, bailing out\n")
		return 1

	# create cumulative set of keys
	all_keys = wig1_index.keys() + wig2_index.keys()
	all_keys = list(set(all_keys))

	print "track type=bedGraph name=test description=avg_test"
	
	#for key in ["chrY"]:
	for key in sorted(all_keys):
		sys.stderr.write("processing {:s}...\n".format(key))

		# create vector of all bases in this chromosome for both files
		fp1 = open("temp1.bin","w+b")
		fp2 = open("temp2.bin","w+b")
		feature_coverage(fp1,key,wig1,wig1_index)
		feature_coverage(fp2,key,wig2,wig2_index)
		fp1.seek(0)
		fp2.seek(0)

		# read through files simultaneously. build a vector of their means
		sys.stderr.write("computing base means...\n")
		chrom_mean = base_means(fp1,fp2,wig1_coverage,wig2_coverage)
		
		fp1.close()
		fp2.close()

		sys.stderr.write("smoothing...\n")
		# chrom_mean = smooth_vector(chrom_mean,20)
		
		# find mean of the two for each base. round each base up so we don't lose anything.
		#sys.stderr.write("mean...\n")
		#chrom_mean = fold_vector(wig1_chrom,wig2_chrom)

		# now we have to bin this down to make a new wig format for this thing
		print_bedgraph(key,chrom_mean)


	#/
	
	os.remove("temp1.bin")
	os.remove("temp2.bin")

#/

def print_bedgraph(chrom,vector):
	num_points = len(vector)
	# out_format = "{:s}\t{:d}\t{:d}\t{:d}\n"
	depth = 0
	depth_last = 0
	start = 0
	i = 0

	depth = vector[i]
	for i in range(1,num_points):
		if depth != depth_last:
			# 
			# found new depth, print out last one
			sys.stdout.write("\t".join(map(str,[chrom,start,i-1,depth_last])) + "\n")
			start = i-1

		depth_last = depth
		depth = vector[i]

	sys.stdout.write("\t".join(map(str,[chrom,start,i,depth_last])) + "\n")

def smooth_vector(v,w):
	num_points = len(v)
	half_w = w/2
	sum = 0
	mean_queue = []
	vout = []
	
	print num_points
	
	for i in range(num_points):
		
		if len(mean_queue) >= half_w:
			vout.append(np.array(mean_queue).mean())
		
		if len(mean_queue) == w:
			mean_queue.pop(0)
		
		mean_queue.append(v[i])
		
		if i % 100000 == 0:
			print len(vout)
	
	while len(mean_queue) > half_w:
		mean_queue.pop(0)
		vout.append(np.array(mean_queue).mean())
		
	return vout

def base_means(fp1,fp2,norm1,norm2):
	vout = []
	v1 = 0
	v2 = 0
	mean = 0
	
	fp1.seek(0,2)
	length_1 = fp1.tell()/2
	fp2.seek(0,2)
	length_2 = fp2.tell()/2
		
	fp1.seek(0)
	fp2.seek(0)
	
	num_points = max(length_1,length_2)
	min_points = min(length_1,length_2)
	temp = 0

	# continue...
	for i in range(min_points):
		v1 = unpack("<H",fp1.read(2))[0] / norm1
		v2 = unpack("<H",fp2.read(2))[0] / norm2
		# round to 6 decimal places
		mean = (v1+v2)/2.0
		mean = round(mean*1000000)/1000000.0
		
		vout.append(mean)

	if min_points != num_points:
		# which one is longer?
		if length_1 > length_2:
			fp_long = fp1
			norm_long = norm1
		else:
			fp_long = fp2
			norm_long = norm2

		for i in range(min_points,num_points):
			v1 = unpack("<H",fp_long.read(2))[0] / norm_long
			v1 = round(v1*1000000)/1000000.0
			vout.append(v1)

	return vout

def fold_vector(v1,v2):
	vout = []
	num_points = max(len(v1),len(v2))
	min_points = min(len(v1),len(v2))
	temp = 0


	# continue...
	for i in range(min_points):
		if v1[i] != 0 and v2[i] != 0:
			temp = log(v2[i]/(v1[i]*1.0),2)
		else:
			temp = 0

		temp  = ceil(temp*1000)
		vout.append(temp/1000.0)

	if min_points != num_points:
		for i in range(min_points,num_points):
			vout.append(0)

	return vout

	
#
# feature_coverage
#
# this function pushes an entire chromosome, at 1 base resolution, into
# a binary file open for writing, fpb.
#
def feature_coverage(fpb,chrom,wig_file,index):

	depth = 0

	# open wig file
	try:
		fp = open(wig_file,"r")
	except IOError,e:
		return 1	

	# skip to chrom
	fp.seek(index[chrom][1])

	# read until end of file or until a new chrom comes up
	line = fp.readline()
	while(line):
		lline = line.strip().split("\t")
		
		if lline[0] == chrom:
			base_first = int(lline[1])
			base_last = int(lline[2])+1
			if base_first > 0:
				base_first += 1
			
			depth = int(lline[3])
			if depth > 0:
				depth -= NOISE_SHIFT
			
			# write each base's coverage out as unsigned 16bit ints
			for i in range(base_first,base_last):
				fpb.write(pack("<H",depth))
			
		line = fp.readline()

	fp.close()

	return 0

	
#
# index_wig
#
# this function indexes the bedgraph file, file, by identifying the file position
# in bytes of the first row for each feature (chromosome) as well as the length
# of the feature (the second value in the last row of a feature's coverage data).
# The total coverage in a file is summed and returned along with the index.
#
def index_wig(file):
	index = {}
	line_last = ""
	chrom = ""
	chrom_start = 0
	coverage = 0

	try:
		fp = open(file,"r")
	except IOError,e:
		return None

	line = ""

	# skip header
	line = fp.readline()
	fpos_last = fp.tell()
	line = fp.readline()

	while(line):
		ll = line.strip().split("\t")
		coverage += int(ll[3])
		if chrom != ll[0]:
			if len(chrom) > 0:
				#
				# store what we learned about this chromosome
				#
				index[chrom] = [int(line_last[2])+1,chrom_start]
			#/

			chrom = ll[0]
			chrom_start = fpos_last
		#/

		line_last = [ll[i] for i in range(len(ll))]
		fpos_last = fp.tell()
		line = fp.readline()
	#/

	# don't forget about the last one!
	index[chrom] = [int(line_last[2])+1,chrom_start]

	fp.close()

	return index,coverage
#/

#
# round
#
# Round a float up or down
#
def round(x):
	x_rem = x - floor(x)
	
	if x_rem >= 0.5:
		return ceil(x)
	
	return floor(x)


if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))


