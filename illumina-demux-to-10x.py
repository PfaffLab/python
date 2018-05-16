#!/usr/bin/python
#==============================================================================
# prep-bd-reads.py
#
# Shawn Driscoll
# 20170807
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script converts the illumina demux files to the 10x format 
# expected in the pachter lab pipeline
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser, isdir
from collections import defaultdict
from time import localtime
import subprocess as sp
from os import system, listdir, mkdir

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	sample_idx = ["ATCGCTCC","CCGTACAG","GATAGGTA","TGACTAGT"]
	index_format = "read-I1_si-{}_lane-001-chunk-001.fastq"
	read_format = "read-RA_si-{}_lane-001-chunk-001.fastq"

	# get list of files in the path.
	flist = listdir(args.fastq_path)

	# get list of sample names and make set of 'R1' files
	samples = {}

	for f in flist:
		tmp = f.split("_")
		sid = []
		for i in range(len(tmp)):
			if re.search("S[0-9]", tmp[i]):
				break
			sid.append(tmp[i])

		tmp = "_".join(sid)
		sid = tmp

		if sid not in samples:
			samples[sid] = []

		if re.search("\_R1\_", f):
			samples[sid].append(f)

	# 
	# now work through each sample
	#

	# 
	# make sure each sample only has 4 files
	#
	hit_err = False
	for sid in samples.keys():
		if len(samples[sid]) > 4:
			error_message("Sample {} has more than 4 split files. Help me?")
			hit_err = True

	if hit_err:
		return 1


	for sid in samples.keys():
		# create output folder
		sid_dir = "{}_fastq".format(sid)
		if not isdir(sid_dir):
			try:
				message("Creating folder {}".format(sid_dir))
				mkdir(sid_dir)
			except OSError as e:
				sys.stderr.write("OSError({}): {}".format(e.errno, e.strerror))
				break
		else:
			message("Folder {} exists. Writing into it.".format(sid_dir))

		# get input read file names
		for i in range(len(samples[sid])):
			rleft = samples[sid][i]
			# make right read file name
			rright = re.sub("\_R1\_", "_R2_", rleft)
			translate(args.fastq_path, sid_dir, rleft, rright, index_format.format(sample_idx[i]), read_format.format(sample_idx[i]), args.b, args.u)

	return 0

def translate(fastq_path, out_path, rleft, rright, index_name, reads_name, bc_len, mi_len):

	lcount = 0
	rcount = 0
	left = []
	right = []

	if re.search("\.gz$", rleft):
		# need to gunzip this guy
		p1 = sp.Popen("gunzip -c {}/{}".format(fastq_path, rleft).split(), stdout=sp.PIPE)
		fin1 = p1.stdout
	else:
		fin1 = open("{}/{}".format(fastq_path, rleft), "r")
	
	if re.search("\.gz$", rright):
		# need to gunzip this guy
		p2 = sp.Popen("gunzip -c {}/{}".format(fastq_path, rright).split(), stdout=sp.PIPE)
		fin2 = p2.stdout
	else:
		fin2 = open("{}/{}".format(fastq_path, rright), "r")
	
	message("parsing {} and {}".format(rleft, rright))
	message("writing {} and {}".format(index_name, reads_name))

	with open("{}/{}".format(out_path, index_name), "w") as fout_index, open("{}/{}".format(out_path, reads_name), "w") as fout_reads:
			
		# make loop runs over first mate file
		for szl in fin1:
			lcount += 1
			szl = szl.strip()
			left.append(szl)
			szl2 = (fin2.readline()).strip()
			right.append(szl2)
					
			if lcount==4:
				rcount += 1

				if (rcount % 1e5) == 0:
					progress_message("parsed {} reads".format(rcount))

				# parse out the cell barcode and the umi
				
				left_read = left[1]
				left_qc = left[3]

				well_bc = left_read[0:bc_len]
				well_qc = left_qc[0:bc_len]
				
				mi = left_read[bc_len:(bc_len+mi_len)]
				mi_qc = left_qc[bc_len:(bc_len+mi_len)]

				#
				# write the cellbarcode to the index fastq
				#
				fout_index.write(left[0]+"\n")
				fout_index.write(well_bc+"\n")
				fout_index.write("+\n")
				fout_index.write(well_qc+"\n")

				#
				# write the read to the reads file
				#
				fout_reads.write(right[0]+"\n")
				fout_reads.write(right[1]+"\n")
				fout_reads.write(right[2]+"\n")
				fout_reads.write(right[3]+"\n")

				#
				# write the UMI
				#
				fout_reads.write(right[0]+"\n")
				fout_reads.write(mi+"\n")
				fout_reads.write("+\n")
				fout_reads.write(mi_qc+"\n")

				lcount = 0
				right = []
				left = []

	fin1.close()
	fin2.close()
	
	progress_message("parsed {} reads".format(rcount))	
	sys.stderr.write("\n")

	cmd = "gzip -f {}/{}".format(out_path, index_name)
	system(cmd)
	cmd = "gzip -f {}/{}".format(out_path, reads_name)
	system(cmd)


def main2(args):

	# variables
	left = []
	right = []
	lcount = 0
	rcount = 0

	sample_idx = ["ATCGCTCC","CCGTACAG","GATAGGTA","TGACTAGT"]

	index_format = "read-I1_si-ATCGCTCC_lane-001-chunk-001.fastq"
	read_format = "read-RA_si-ATCGCTCC_lane-001-chunk-001.fastq"
	
	if not isfile(args.mate1):
		error_message("Input file {} does not exist".format(args.mate1))
	
	if not isfile(args.mate2):
		error_message("Input file {} does not exist".format(args.mate2))
	
	#
	# read both files and handle each read one at a time
	#

	if re.search("\.gz$", args.mate1):
		# need to gunzip this guy
		p1 = sp.Popen("gunzip -c {}".format(args.mate1).split(), stdout=sp.PIPE)
		fin1 = p1.stdout
	else:
		fin1 = open(args.mate1, "r")
	
	if re.search("\.gz$", args.mate2):
		# need to gunzip this guy
		p2 = sp.Popen("gunzip -c {}".format(args.mate2).split(), stdout=sp.PIPE)
		fin2 = p2.stdout
	else:
		fin2 = open(args.mate2, "r")
	

	with open(index_format, "w") as fout_index, open(read_format, "w") as fout_reads:
			
		# make loop runs over first mate file
		for szl in fin1:
			lcount += 1
			szl = szl.strip()
			left.append(szl)
			szl2 = (fin2.readline()).strip()
			right.append(szl2)
					
			if lcount==4:
				rcount += 1

				if (rcount % 1e5) == 0:
					progress_message("parsed {} reads".format(rcount))

				# parse out the cell barcode and the umi
				
				left_read = left[1]
				left_qc = left[3]

				well_bc = left_read[0:8]
				well_qc = left_qc[0:8]
				
				mi = left_read[8:16]
				mi_qc = left_qc[8:16]

				#
				# write the cellbarcode to the index fastq
				#
				fout_index.write(left[0]+"\n")
				fout_index.write(well_bc+"\n")
				fout_index.write("+\n")
				fout_index.write(well_qc+"\n")

				#
				# write the read to the reads file
				#
				fout_reads.write(right[0]+"\n")
				fout_reads.write(right[1]+"\n")
				fout_reads.write(right[2]+"\n")
				fout_reads.write(right[3]+"\n")

				#
				# write the UMI
				#
				fout_reads.write(right[0]+"\n")
				fout_reads.write(mi+"\n")
				fout_reads.write("+\n")
				fout_reads.write(mi_qc+"\n")

				lcount = 0
				right = []
				left = []

	fin1.close()
	fin2.close()
	
	progress_message("parsed {} reads".format(rcount))	
	sys.stderr.write("\n")

	cmd = "gzip -f {}".format(index_format)
	os.system(cmd)
	cmd = "gzip -f {}".format(read_format)
	os.system(cmd)

	return 0

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz, show_time=True):
	if show_time:
		sys.stderr.write("[{}] {}\n".format(time_string(), sz))
	else:
		sys.stderr.write("{}\n".format(sz))

def progress_message(sz):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r[{}] {}".format(time_string(), sz))

def message_mp(sz, name, lock):
	lock.acquire()
	sys.stderr.write("[{}] {}\n".format(name, sz))
	lock.release()
	return


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Converte normal Illumina demultiplexed 10x reads into the index/read format that cellrange exports.")
parser.add_argument('fastq_path', type=str, help="Path to fastq files")

parser.add_argument('-b', type=int, default=16, action="store", 
	help="Cell barcode length. 14 or 16 [16]")
parser.add_argument('-u', type=int, default=10, action="store", 
	help="UMI barcode length. [10]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

