#!/usr/bin/env python
#==============================================================================
# short-read-sim.py
#
# Shawn Driscoll
# 20120824
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Reads in a FASTA file and generates short reads (illumina fastq style)
# for each sequence in a random fashion. Reads can be single or paired-end.
#
# Updates:
# 20130403
# - changed name of the script to short-read-sim.py 
# - adding in an optional error rate to create base call errors
#==============================================================================

import sys, argparse, math, re
from subprocess import Popen
from random import gauss, random, sample
from scipy.stats import norm
import numpy as np
import numpy.random as npr

import rpy2.robjects as robjects
r = robjects.r

#==============================================================================
# main
#==============================================================================
def main(args):

	# variables
	fp = None
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

	fai = args.infile + ".fai"
	dfai = {}

	# 
	# check input file
	#
	if not file_exists(args.infile):
		sys.stderr.write("Error: input file doesn't exist!\n")
		return(1)
	#}

	#
	# check for fasta index
	#
	if not file_exists(fai):
		sys.stderr.write("> Cannot find fasta index, rebuilding now...\n")
		p1 = Popen("samtools faidx {:s}".format(args.infile).split())
		p1.wait()
	#}

	# hash fasta index
	sys.stderr.write("> Hashing fasta index...")
	fin = open(fai, "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		dfai[ll[0]] = list(map(int, ll[1:]))
	#}
	fin.close()
	sys.stderr.write("done\n")

	#
	# set minimum reference length based on read length or fragment length
	if args.make_paired:
		min_seq_len = args.frag_len * 3
	else:
		min_seq_len = args.read_length * 3
	#}

	#
	# check if there is a map file
	#
	if args.map is not None or args.make_random_map > 0:
		
		if args.map is not None:
			# hash map file
			using_map = True
			fin = open(args.map, "r")
			for szl in fin:
				ll = szl.strip().split("\t")
				if len(ll) > 1:
					dmap[ll[0]] = float(ll[1])
				else:
					dmap[ll[0]] = args.coverage
				#}
			#}
	
			fin.close()
		else:

			#
			# Generate a random set of N features where N = args.make_random_map.
			#
			
			# make a new dict of fasta ref ids that are longer than the minimum
			# sequence length			
			dfai_temp = {}
			for rid in dfai.keys():
				if dfai[rid][0] > min_seq_len:
					# keep it
					dfai_temp[rid] = dfai[rid][0]
				#}
			#}
			
			if len(dfai_temp.keys()) < args.make_random_map:
				# not enough references are long enough to make the requested number of random
				# refs
				sys.stderr.write("Warning: only {:d} references are long enough to sample\n")
				args.make_random_map = len(dfai_temp.keys())
			#}
			
			# make a random set
			rset = sample(dfai_temp.keys(), args.make_random_map)
			
			dcounts = []
			if args.pseudo_map:
				dcounts = make_multinom_counts(args.total_reads, args.make_random_map)
			
			
				
			dmap = {}
			i = 0
			for rid in rset:
				if args.prob_map:
					dmap[rid] = dfai[rid][0]
					
				elif args.pseudo_map:
					if args.make_paired:
						dmap[rid] = dcounts[i]*args.read_length*2.0/dfai[rid][0]
					else:
						dmap[rid] = dcounts[i]*args.read_length*1.0/dfai[rid][0]
					i += 1
						
				else:
					dmap[rid] = args.coverage
				#}
			#}

		#}
		
		# check list against the fasta index
		for rid in dmap.keys():
			if rid in dfai:
				# check length against minimum length
				if dfai[rid][0] < min_seq_len:
					# remove feature
					dmap.pop(rid, None)
					sys.stderr.write("Warning: removing {:s} from simulation because it's too short ({:d})\n".format(rid, dfai[rid][0]))
				#}
			#}
		#}

		if args.prob_map:
			# values dmap are intended to be proportions of reads. first we have to 
			# make them sum to one and convert them to total read counts
			n_sum = 0
			for rid in dmap.keys():
				n_sum += dmap[rid]
			n_ratio = 1.0/n_sum

			for rid in dmap.keys():
				dmap[rid] = dmap[rid]*n_ratio * args.total_reads
				if dmap[rid] < 1:
					dmap[rid] = 1

				# scale by read length (or read length x2 for paired)
				if args.make_paired:
					dmap[rid] *= args.read_length*2
				else:
					dmap[rid] *= args.read_length

				# now the total base count is converted to coverage based on the
				# feature's length (from the fasta index)
				dmap[rid] = dmap[rid]/dfai[rid][0]
			#}
		#}

	elif args.make_random_map == -1:
		# for this mode all transcripts that are long enough to be sampled will be used and a semi-realistic count
		# distribution will be randomly generated.
		sys.stderr.write("> generating multinomial count distribution...")
		# first figure out which references can generate reads. this will be the dmap
		dmap = {}
		count_prob = []
		for rid in sorted(dfai.keys()):
			if dfai[rid][0] > min_seq_len:
				# keep it
				dmap[rid] = dfai[rid][0]
				count_prob.append(10**r.rnorm(1, 2.5, 1)[0]-0.01)
				if count_prob[-1] < 0:
					count_prob[-1] = 0
		
		# generate random multinomial count distribution
		count_levels = r.rmultinom(1, args.total_reads, count_prob)
		# copy count_levels into the dmap and convert to coverages (they get converted back into counts
		# in the read generating functions)
		i =0
		for rid in sorted(dmap.keys()):
			if args.make_paired:
				dmap[rid] = count_levels[i]*args.read_length*2.0/dfai[rid][0]
			else:
				dmap[rid] = count_levels[i]*args.read_length*1.0/dfai[rid][0]
			i += 1 
		
		sys.stderr.write("done\n")

	else:
		# no map file so we'll take all features
		for rid in dfai.keys():
			dmap[rid] = args.coverage
		#}
	#}
	
	#
	# split off extention from input file
#	ltemp = os.path.basename(args.infile).split(".")
#	fout_stub = ".".join(ltemp[0:(len(ltemp)-1)])

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
	#
	fpin = open(args.infile,"r")

	# loop through reference id's to be sampled
	for rid in dmap.keys():
		# fetch sequence from FASTA for current reference
		if rid in dfai:

			if dfai[rid][0] >= min_seq_len and dmap[rid] > 0:

				sequence = fetch_ref_seq(fpin, dfai, rid).upper()

				# 
				# initalize lists
				lreads = []
				rreads = []
				if args.make_paired:
					lheaders = [[], []]
				else:
					lheaders = []
				header = rid
				coverage = dmap[rid]

				# 
				# make reads from current sequence
				if args.make_paired:
					result = make_paired_reads(lreads, rreads, lheaders, header, sequence, args.read_length, coverage, args.frag_len, args.frag_sdev, args.error_rate)
				else:
					result = make_reads(lreads, lheaders, header, sequence, args.read_length, coverage, args.error_rate)

				# check return value
				if result == 2:
					# show warning about skipping a sequence due to length
					sys.stderr.write("warning: skipping {:s} because it's too short ({:d})\n".format(header, len(sequence)))
				elif result != 0:
					# something wrong happened
					fpin.close()
					if args.make_paired:
						fp1.close()
						fp2.close()
					else:
						fp.close()
					sys.exit(1)

				# write reads out to file
				if result == 0:
					if args.make_paired:
						print_reads(fp1, lreads, lheaders[0])
						print_reads(fp2, rreads, lheaders[1])
					else:
						print_reads(fp, lreads, lheaders)
			else:
				if dmap[rid] == 0:
					pass
				else:
					sys.stderr.write("warning: skipping {:s} because it's too short ({:d})\n".format(rid, dfai[rid][0]))
		else:
			sys.stderr.write("warning: reference is not in the fasta index ({:s})\n".format(rid))
		#}
	#}

	fpin.close()

	if args.make_paired:
		fp1.close()
		fp2.close()
	else:
		fp.close()


	return 0

#==============================================================================
# make_reads
# Creates random reads from single sequence and passes them to a list
#==============================================================================
def make_reads(lout, lheaders, seq_name, seq, read_length, cov, err):

	num_reads = 0
	made_reads = 0
	ppos = 0
	made_reads = 0
	flipped = False
	seq_mean = 0
	seq_sd = 0
	rheader = ""
	bases = ['A', 'C', 'T', 'G']
	prob = 0
	e_base = ""
	e_string = ""
	read = ""
	read_new = ""
	seq_rev = reverse_comp(seq)

	# do nothing if the sequence length less than read_length
	if len(seq) < (read_length*2):
		return 2

	# get sequence length
	seqlen = len(seq)
	# center (mean) of sequence
	seq_mean = (seqlen-read_length)/2
	# set sd for random read selection about the center of the sequence such that
	# reads from the end of the sequence are within 1SD
	seq_sd = seq_mean
	temp = ""
	
	if cov > 0:
		# calculate number of reads to create
		num_reads = int(math.ceil((seqlen * cov)/(read_length*1.0)))
	else:
		num_reads = 1
	#}

	# determin maximum start position within sequence so that any random read 
	# doesn't go past the end of the sequence
	max_pos = seqlen - read_length - 1
	min_pos = 0


	# main loop
	while made_reads < num_reads:
		# pull reads from the sequence in sort of a normally distributed fashion biasing
		# towards the center of the sequence
		ppos = gauss(seq_mean, seq_sd)
		if ppos > (max_pos+read_length/2) or ppos < -(read_length/2):
			pass
		else:
			if ppos < 0:
				ppos = 0
			elif ppos > max_pos:
				ppos = max_pos

			e_string = ""
			ppos = int(math.floor(ppos))

			# 50/50 chance of read being reversed
			# flipped = norm.cdf(gauss(0, 1)) > 0.5
			flipped = r.runif(1, 0, 1)[0] > 0.5

			# get the read
#			if flipped:
#				read = seq_rev[ppos:(ppos+read_length)]
				# need to revise ppos
#				temp = seqlen - ppos - read_length
#				ppos = temp 
#			else:
			read = seq[ppos:(ppos+read_length)]

			if len(read) == read_length:

				if err > 0:
					# generate base errors
					read_new, e_string = generate_base_errors(read, err)
					read = read_new
					
				if flipped:
					read = reverse_comp(read)

				lout.append(read)

				# make the header
				if len(e_string) == 0:
					e_string = "NE"

				rheader = "@{:s}_{:d}_{:d}_{:d}_{:d}_{:s}".format(seq_name, ppos+1, ppos+read_length, int(flipped), made_reads, e_string)

				lheaders.append(rheader)
				made_reads += 1


	return 0

#==============================================================================
# make_reads
# Creates random reads from single sequence and passes them to a list
#==============================================================================
def make_paired_reads(lleft, lright, lheaders, seq_name, seq, read_length, cov, frag_len, frag_sdev, err):

	num_reads = 0
	made_reads = 0
	ppos = 0
	ppos2 = 0
	made_reads = 0
	flipped = False
	seq_mean = 0
	seq_sd = 0
	rheader = ""
	bases = ['A', 'C', 'T', 'G']
	prob = 0
	e_base = ""
	e_string1 = ""
	e_string2 = ""
	read1 = ""
	read2 = ""
	read_new = ""
	seq_rev = reverse_comp(seq)
	flen = 0

	# do nothing if the sequence length less than read_length
	if len(seq) < (frag_len*2):
		return 2

	# get sequence length
	seqlen = len(seq)
	# center (mean) of sequence
	seq_mean = (seqlen-frag_len)/2
	# set sd for random read selection about the center of the sequence such that
	# reads from the end of the sequence are within 1SD
	seq_sd = seq_mean
	temp = ""
	
	if cov > 0:
		# calculate number of reads to create
		num_reads = int(math.ceil((seqlen * cov)/(read_length*2.0)))
	else:
		num_reads = 1
	#}


	# main loop
	while made_reads < num_reads:
		# calc fragment length
		
		flen = 0
		while flen < (read_length+1):
			# make sure fragment length is never less than the read length+1
			flen = int(math.floor(gauss(frag_len, frag_sdev)))

		max_pos = seqlen-flen-1

		# pull reads from the sequence in sort of a normally distributed fashion biasing
		# towards the center of the sequence
		seq_mean = (seqlen-frag_len)*1.0/2
		ppos = gauss(seq_mean, seq_mean)

		if ppos > (max_pos+read_length/2) or ppos < 0-(read_length/2):
			pass
		else:
			if ppos < 0:
				ppos = 0
			elif ppos > max_pos:
				ppos = max_pos

			e_string1 = ""
			e_string2 = ""
			ppos = int(math.floor(ppos))
			ppos2 = (ppos + flen) - read_length

			# 50/50 chance of read being reversed
#			flipped = norm.cdf(gauss(0, 1)) > 0.5
			flipped = r.runif(1, 0, 1)[0] > 0.5

			# pull the fragment
			frag = seq[ppos:(ppos+flen)]
			read1 = frag[0:read_length]
			read2 = reverse_comp(frag[(flen-read_length):])

			# get the reads
#			if flipped:
#				read1 = seq_rev[ppos:(ppos+read_length)]
#				read2 = reverse_comp(seq_rev[ppos2:(ppos2+read_length)])
#			else:
#				read1 = seq[ppos:(ppos+read_length)]
#				read2 = reverse_comp(seq[ppos2:(ppos2+read_length)])

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
					

				lleft.append(read1)
				lright.append(read2)

				# make the header
				if len(e_string1) == 0:
					e_string1 = "NE"
				if len(e_string2) == 0:
					e_string2 = "NE"

				rheader = "@{:s}_{:d}_{:d}_{:d}|{:d}_{:d}_{:s}_{:s}".format(seq_name, ppos+1, ppos+flen, int(flipped), int(not flipped), made_reads, e_string1, e_string2)
				
				lheaders[0].append(rheader + "/1")
				lheaders[1].append(rheader + "/2")

				made_reads += 1

	return 0

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
	seq = re.sub('\n', '', fin.read(read_length))

	return seq.strip()


def make_multinom_counts(depth, num_features):
	
	count_prob = [0 for i in range(num_features)]
	for i in range(num_features):
		count_prob[i] = 10**r.rnorm(1, 2.5, 1)[0]-0.01
		
	# generate random multinomial count distribution
	counts = r.rmultinom(1, args.total_reads, count_prob)
	
	return list(counts)	


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
parser.add_argument('-c',dest="coverage",default=2,type=float, help="Coverage per sequence. (default: 2)")
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
parser.add_argument('-N', dest="total_reads", action="store", default=1000000, type=int, 
					help="with -m and --probability-map this value is used as a cap for total read depth (default: 1000000)")
parser.add_argument('--min-read-depth', dest="min_reads", action="store", type=int, default=10, 
					help="with -m and --probability-map this value is used as the minimum read depth (default: 10)")
parser.add_argument('--random-map', dest="make_random_map", type=int, action="store", default=0, 
					help="Generate a random map of references to sample from the FASTA index of references. If you enter -1 the full list of transcripts is used however a random, semi-realistic count distribution will be generated. Some transcripts my contain zero counts.")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
