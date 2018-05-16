#!/usr/bin/env python
#==============================================================================
# short-read-sim-mp.py
#
# Shawn Driscoll
# 20120824
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Reads in a FASTA file and generates short reads (illumina fastq style)
# for each sequence in a random fashion. Reads can be single or paired-end.
# This version can run several parallel threads for maximum quickness.
#
# Updates:
# 20130403
# - changed name of the script to short-read-sim.py 
# - adding in an optional error rate to create base call errors
#==============================================================================

import sys, argparse, math
from subprocess import Popen
from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
import numpy.random as npr
import multiprocessing as mp

# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================
def main(args):

	# variables
	fp1 = None
	fp2 = None
	sequence = ""
	cline = ""
	lreads = []
	rreads = []
	lheaders = []
	header = ""
	result = 0
	seq_temp = ""
	fout_stub = args.stub
	ltemp = []
	dmap = {}
	using_map = False
	min_seq_len = 0
	made_reads = 0
	num_keys = 0
	num_threads = 0
	num_tasks = 0
	total_seqs = 0
	i = 0
	j = 0
	
	result = []
	
	# queues and lock for multithread workers
	qtasks = mp.Queue()
	qresult = mp.Queue()
	lock = mp.Lock()	

	# make fasta index file name from input fasta name
	fai = args.infile + ".fai"
	dfai = {}

	# make sure the user didn't specify mutually exclusive options
	if args.pseudo_map and args.prob_map:
		sys.stderr.write("Error: you can't specify both --probability-map and --pseudo-map\n")
		return 1

	# 
	# check input file
	if not file_exists(args.infile):
		sys.stderr.write("Error: input file doesn't exist!\n")
		return(1)

	#
	# check for fasta index
	if not file_exists(fai):
		sys.stderr.write("> cannot find fasta index, rebuilding now...\n")
		p1 = Popen("samtools faidx {:s}".format(args.infile).split())
		p1.wait()

	#
	# set minimum reference length based on read length or fragment length
	if args.make_paired:
		min_seq_len = args.frag_len * 3
	else:
		min_seq_len = args.read_length * 3
		

	# hash fasta index
	fin = open(fai, "r")
	dfai = {}
	for szl in fin:
		ll = szl.strip().split("\t")
		# only keep sequences that are long enough for this simulation
		if int(ll[1]) >= min_seq_len:
			total_seqs += 1
			dfai[ll[0]] = list(map(int, ll[1:]))

	fin.close()
	sys.stderr.write("> hashed %d sequence ids from fasta index\n" % total_seqs)

	# build the dmap dict which will contain the reference id's to process and their read counts/coverages
	dmap = {}
	
	# if there is a map file we can read it in
	if args.map is not None:
		# hash map file but drop sequences that aren't in dfai
		using_map = True
		fin = open(args.map, "r")
		for szl in fin:
			ll = szl.strip().split("\t")
			if ll[0] in dfai:
				
				if len(ll) > 1:
					# map file contains some information - copy that in
					dmap[ll[0]] = float(ll[1])
				else:
					# default to use the coverage parameter. convert that to 
					# read counts
					dmap[ll[0]] = round(args.coverage*dfai[ll[0]][0]/args.read_length)
					
			else:
				sys.stderr.write("> warning: dropping %s from map file because it is too short\n" % ll[0])
				
		fin.close()
	
	elif args.make_random_map > 0:
		# select random reference ids
		if total_seqs < args.make_random_map:
			# not enough references are long enough to make the requested number of random
			sys.stderr.write("> warning: only {:d} references are long enough to sample\n".format(total_seqs))
			args.make_random_map = total_seqs
		
		# make a random set
		rset = sample(dfai.keys(), args.make_random_map)
		
		for rid in rset:
			dmap[rid] = round(args.coverage*dfai[rid][0]/args.read_length)
	
	else:
		# draw from all references		
		for rid in dfai.keys():
			dmap[rid] = round(args.coverage*dfai[rid][0]/args.read_length)
	
	##
	## now deal with what type of count distribution we want
	##
	
	if args.pseudo_map:
		dcounts = make_multinom_pseudo_counts(args.total_reads, len(dmap.keys()))
		i = 0
		for rid in dmap.keys():
			if args.make_paired:
				dmap[rid] = dcounts[i]
			else:
				dmap[rid] = dcounts[i]
			i += 1
				
	elif args.prob_map:
		if using_map:
			# use values pulled from map file as a probability distribution and then 
			# generate a random multinomial distribtuion of counts from it
			probs = [dmap[rid] for rid in dmap.keys()]
			dcounts = make_multinom_counts(args.total_reads, probs)
		else:
			dcounts = make_multinom_uniform_counts(args.total_reads, len(dmap.keys()))
				
		i = 0
		for rid in dmap.keys():
			if args.make_paired:
				dmap[rid] = dcounts[i]
			else:
				dmap[rid] = dcounts[i]
			i += 1
	
	elif args.length_prob:
		probs = [dfai[rid][0] for rid in dmap.keys()]

		dcounts = make_multinom_counts(args.total_reads, probs)
		
		i = 0
		for rid in dmap.keys():
			if args.make_paired:
				dmap[rid] = dcounts[i]
			else:
				dmap[rid] = dcounts[i]
			i += 1		
	
	# else we just use what's in the dmap already
	
	if args.save_map:
		fout = open("srs-map.txt", "w")
		for rid in sorted(dmap.keys()):
			fout.write("%s\t%d\n" % (rid, dmap[rid]))
		fout.close()
	
	#
	# open output file(s)
	if args.make_paired:
		fp1 = open(fout_stub + ".1.fq","w")
		fp2 = open(fout_stub + ".2.fq","w")
	else:
		fp = open(fout_stub + ".fq","w")

	sys.stderr.write("> Generating reads...\n")

	#
	# open fasta file and get to work
	fpin = open(args.infile,"r")

	num_tasks = 0

	# make workers
	for i in range(args.num_threads):
		mp.Process(target=worker, args=(qtasks, qresult, lock)).start()
		
	# load up the tasks in the queue so the workers have something to d
	
	klist = dmap.keys()
	num_k = len(klist)
	i = 0
	while i < num_k:
		
		# send up to 100 tasks at a time so the buffer doesn't get too full
		num_tasks = 0
		for j in range(100):
			rid = klist[i]
					
			if dmap[rid] > 0:
				# pull sequence
				sequence = fetch_ref_seq(fpin, dfai, rid).upper()
				
				# push into queue
				if args.make_paired:
					qtasks.put((make_paired_reads, (rid, sequence, args.read_length, dmap[rid], args.frag_len, args.frag_sdev, args.error_rate)))
				else:
					qtasks.put((make_reads, (rid, sequence, args.read_length, dmap[rid], args.error_rate)))
				
				# count tasks
				num_tasks += 1
				#print "%d: sent %s" % (num_tasks, rid)
			
			i += 1
			if i >= num_k:
				break

#		sys.stderr.write("> processing %d sequences...\n" % num_tasks)
		
		# retrieve results from this batch of tasks and print out to file(s)
		for j in range(num_tasks):
			result = qresult.get()
#			print j, len(result)
			if args.make_paired:
				write_paired_reads(fp1, fp2, result)
			else:
				write_reads(fp, result)
	
	# insert the STOP order at the end of the queue since everything has been sent
	for i in range(args.num_threads):
		qtasks.put('STOP')
	
	fpin.close()

	if args.make_paired:
		fp1.close()
		fp2.close()
	else:
		fp.close()

	return 0

#===============================================================================
# worker
# Working function for each parallel process. This function simply iterates
# through available tasks in qin and dumps output to qout. 
#===============================================================================
def worker(qin, qout, lock):
	for func, args in iter(qin.get, 'STOP'):
#		print "worker"
		result = func(*args)
#		print "got result"
#		lock.acquire()
		qout.put(result)
#		print "put result"
#		lock.release()


#==============================================================================
# make_reads
# Creates random reads from single sequence and passes them to a list
#==============================================================================
def make_reads(seq_name, seq, read_length, num_reads, err):

	made_reads = 0
	ppos = 0
	flipped = False
	rheader = ""
	e_string = ""
	read = ""
	read_new = ""
	
	# reads list, each entry will be a header/read pair in a list
	lreads = []

#	sys.stderr.write("[make_reads] %s %d tick..." % (seq_name, len(seq)))

	# get sequence length
	seqlen = len(seq)
	
	if seqlen < read_length:
		sys.stderr.write("> warning: %s snuck through with length %d\n" % (seq_name, seqlen))
		return []
	
	if num_reads < 1:
#		sys.stderr.write("BAIL\n")
		return []

	# find maximum start position within sequence so that any random read 
	# doesn't go past the end of the sequence
	max_pos = seqlen - read_length - 1
	min_pos = 0

	# main loop
	while made_reads < num_reads:
		# pull reads from the sequence in sort of a normally distributed fashion biasing
		# towards the center of the sequence
		
		ppos = int(round(npr.rand()*max_pos))
		
		if ppos > max_pos:
			ppos = max_pos
		
		e_string = ""

		# 50/50 chance of read being reversed
		flipped = npr.rand() > 0.5

		read = seq[ppos:(ppos+read_length)]

		if len(read) == read_length:

			if err > 0:
				# generate base errors
				read_new, e_string = generate_base_errors(read, err)
				read = read_new
				
			if flipped:
				read = reverse_comp(read)

			# make the header
			if len(e_string) == 0:
				e_string = "NE"

			rheader = "@{:s}_{:d}_{:d}_{:d}_{:d}_{:s}".format(seq_name, ppos+1, ppos+read_length, int(flipped), made_reads, e_string)

			lreads.append([rheader, read])
			
			# increment read count
			made_reads += 1
	
#	print "make_reads made %d reads from %s" % (len(lreads), seq_name)
#	sys.stderr.write("tock\n")
	return lreads

#==============================================================================
# make_reads
# Creates random reads from single sequence and passes them to a list
#==============================================================================
def make_paired_reads(seq_name, seq, read_length, num_reads, frag_len, frag_sdev, err):

	made_reads = 0
	ppos = 0
	flipped = False
	rheader = ""
	e_string1 = ""
	e_string2 = ""
	read1 = ""
	read2 = ""
	read_new = ""
	flen = 0
	
	lreads = []

	# get sequence length
	seqlen = len(seq)
	
	if num_reads < 1:
		return []

	# main loop
	while made_reads < num_reads:
		
		# calc fragment length for current fragment
		flen = 0
		while flen < (read_length+1):
			# make sure fragment length is never less than the read length+1
			flen = int(math.floor(npr.normal(frag_len, frag_sdev)))

		max_pos = seqlen-flen-1

		# ppos = int(round(r.runif(1, 0, 1)[0]*max_pos))
		ppos = int(round(npr.rand()*max_pos))

		if ppos > max_pos:
			ppos = max_pos

		e_string1 = ""
		e_string2 = ""

		# 50/50 chance of read being reversed
		flipped = npr.rand() > 0.5

		# pull the fragment
		frag = seq[ppos:(ppos+flen)]
		read1 = frag[0:read_length]
		read2 = reverse_comp(frag[(flen-read_length):])

		if flipped:
			# for pairs if the read is reversed then the second mate is really the 
			# first mate. no more complimenting of the reads is necessary
			temp = read1
			read1 = read2
			read2 = temp

		if len(read1) == read_length and len(read2) == read_length:

			if err > 0:
				# generate base errors
				read_new, e_string1 = generate_base_errors(read1, err)
				read1 = read_new
				read_new, e_string2 = generate_base_errors(read2, err)
				read2 = read_new
				
			# make the header
			if len(e_string1) == 0:
				e_string1 = "NE"
			if len(e_string2) == 0:
				e_string2 = "NE"

			rheader = "@{:s}_{:d}_{:d}_{:d}_{:d}_{:s}_{:s}".format(seq_name, ppos+1, ppos+flen, int(flipped), made_reads, e_string1, e_string2)
			
			lreads.append([rheader, read1, read2])
			# increment read count
			made_reads += 1

	return lreads

def generate_base_errors(read, err):
	bases = ['A', 'C', 'T', 'G']
	prob = 0
	e_base = ""
	e_string = ""
	read_new = ""
	e_string = ""
	read_length = len(read)

	# x = abs(npr.randn(read_length))
	# y = x * -1
	# prob = 1 - (norm.cdf(x) - norm.cdf(y))
	# prob = list(r.runif(read_length, 0, 1))
#	prob = [random() for i in range(read_length)]
#	prob_ci = np.where(np.array(prob) < err)[0]

	prob = npr.rand(read_length)
	prob_ci = np.where(prob < err)[0]

	if len(prob_ci) == 0:
		return (read, "")

	lread = list(read)
	for i in range(len(prob_ci)):
		e_base = sample(bases, 1)[0]
		while e_base == lread[prob_ci[i]]:
			e_base = sample(bases, 1)[0]

		e_string += "{:d}{:s}>{:s}".format(prob_ci[i]+1, lread[prob_ci[i]], e_base)
		lread[prob_ci[i]] = e_base

	read = "".join(lread)

	return (read, e_string)


def write_paired_reads(fp1, fp2, reads):
	
	n = len(reads)
	read_length = len(reads[0][1])
	qual_string = "".join([chr(30+33) for i in range(read_length)])
	
	for i in range(n):
		fp1.write("{:s}/1\n{:s}\n+\n{:s}\n".format(
					reads[i][0],
					reads[i][1],
					qual_string))
		fp2.write("{:s}/2\n{:s}\n+\n{:s}\n".format(
					reads[i][0],
					reads[i][2],
					qual_string))
	
	return 0

def write_reads(fp, reads):

	n = len(reads)
	read_length = len(reads[0][1])
	qual_string = "".join([chr(30+33) for i in range(read_length)])
	
	for i in range(n):
		fp.write("{:s}\n{:s}\n+\n{:s}\n".format(
					reads[i][0],
					reads[i][1],
					qual_string))
	
#	print "wrote %d reads\n" % n
	
	return 0
	

#==============================================================================
# print_reads
# Prints reads in FASTQ format. All base qualities are set to 40
#==============================================================================
def print_reads(fp, lout, lheaders):

	if len(lout) == 0:
		return 1

	# get read length
	read_length = len(lout[0])
	# create qual string
	qual_string = "".join([chr(30+33) for i in range(read_length)])

	#
	# loop through reads and print out in FASTQ format to file
	for i in range(len(lout)):
		# read header
		fp.write(lheaders[i] + "\n")
		fp.write("{:s}\n".format(lout[i]))
		fp.write("+\n")
		fp.write("{:s}\n".format(qual_string))


	return 0

#
# fetch_ref_seq
# Fetches sequence for specified reference from the open FASTA file, fin
# using the hashed fasta index, fai, as a guide
def fetch_ref_seq(fin, fai, ref_name):

	s_length = fai[ref_name][0]
	s_offset = fai[ref_name][1]

	# compute number of lines the sequence takes up (less a partial line)
	s_lines = s_length / fai[ref_name][2]

	# compute actual in-file length of the sequence. that's the number of lines at 
	# the column 4 length plus however many characters are in the last line
	read_length = s_lines*fai[ref_name][3] + (s_length - s_lines*fai[ref_name][2])

	# read in the sequence and strip out any newlines
	fin.seek(s_offset)
	seq = fin.read(read_length)
	seq = seq.replace("\n", "")
	#seq = re.sub('\n', '', fin.read(read_length))
		
	return seq.strip()


#==============================================================================
# reverse_comp
# Returns the reverse complement of a sequence
#==============================================================================
def reverse_comp(seq):

	# reverse compliment lookup table
	rev_cmp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
	seqlen = len(seq)
	seqout = ""
	i = 0

	if seqlen == 0:
		return ""

	i = seqlen
	while i > 0:
		i -= 1
		seqout += rev_cmp[seq[i]]

	return seqout

def make_multinom_pseudo_counts(depth, num_features):
	
	count_prob = npr.normal(2.5, 1, num_features)
	count_prob = 10**count_prob-0.01
	ci = np.where(count_prob < 0)[0]
	if len(ci) > 0:
		for i in range(len(ci)):
			count_prob[ci[i]] = 0
	
	# generate random multinomial count distribution
	# counts = r.rmultinom(1, depth, list(count_prob))
	counts = npr.multinomial(depth, count_prob/count_prob.sum(), size=1)[0]
	
	return list(counts)	

def make_multinom_uniform_counts(depth, num_features):
	
	# generate random uniform counts
	# count_prob = [random() for i in range(num_features)]
	count_prob = npr.rand(num_features)
		
	# generate random multinomial count distribution
	# counts = r.rmultinom(1, depth, count_prob)
	counts = npr.multinomial(depth, count_prob/count_prob.sum(), size=1)[0]
	
	return list(counts)	

def make_multinom_counts(depth, prob):
	# generate random multinomial count distribution from supplied probability vector
#	counts = r.rmultinom(1, depth, prob)
	prob = np.array(prob, dtype=float)
	counts = npr.multinomial(depth, prob/prob.sum(), size=1)
	return list(counts[0])	
	

def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Short read (SE or PE) simulator for transcriptome references.")
parser.add_argument('infile',type=str,help="FASTA input file")
parser.add_argument('-c',dest="coverage",default=2,type=float, help="Coverage per sequence. Ignored if --probability-map or --pseudo-map are specified. (default: 2)")
parser.add_argument('-p', dest="num_threads", default=1, type=int, help="Number of threads to run. (default: 1)")
parser.add_argument('-o', dest="stub", default="reads", type=str, help="Output reads file stub (default: reads)")
parser.add_argument('-l',dest="read_length",default=50,type=int,help="Read length. (default: 50)")
parser.add_argument("-e", dest="error_rate", default=0.01, type=float, help="Base error rate. Illumina typical could be 0.001 (default: 0.01)")
parser.add_argument('--paired',dest="make_paired",action="store_const",default=False,const=True,help="Generate paired-end reads (default: False)")
parser.add_argument('-fm',dest="frag_len",type=int,action="store",default=200,
				help="(with --paried) Set mean fragment length of pairs (default: 200)")
parser.add_argument('-fs',dest="frag_sdev",type=int,action="store",default=60,
					help="(with --paried) Set fragment length standard deviation (default: 60)")
parser.add_argument('-m', dest="map", type=str, action="store", default=None, 
					help="Provide a map index indicating which references to draw reads from (default: All)")
parser.add_argument('--probability-map', dest='prob_map', action="store_const", default=False, const=True, 
					help="with -m, map contains probability of reads per feature and not coverage rates (default: off)")
parser.add_argument('--pseudo-map', dest='pseudo_map', action="store_const", default=False, const=True, 
					help="with --random-map, generate pseudo-realistic count distribution (default: off)")
parser.add_argument('--length-probs', dest='length_prob', action="store_const", default=False, const=True, 
					help="assign probability of read assignment based on feature length (default: off)")
parser.add_argument('-N', dest="total_reads", action="store", default=1000000, type=int, 
					help="with -m, --probability-map, --pseudo-map or --length-probs this value is used as a cap for total read depth (default: 1000000)")
parser.add_argument('--random-map', dest="make_random_map", type=int, action="store", default=0, 
					help="Generate a random map of references to sample from the FASTA index of references. Some transcripts my contain zero counts depending on read depth (-N).")
parser.add_argument('--save-map', dest="save_map", action="store_const", const=True, default=False,
				help="Save the dynamically generated transcript/coverage map to file as srs-map.txt (default: False)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
