#!/usr/bin/env python

import sys, pysam, argparse
from math import floor
import hashlib
import numpy as np
from numpy.random import rand, normal
from random import sample
import re

def main(args):
	
	flast = 0
	hfai = {}
	offset = 0
	rev = False
	rl = args.read_length
	fl = 0
	rbuff = ""
	rcount = 0
	e_rate = 0.01
	rmm = 0
	
	qname = ""
	quals = "".join([chr(30+33) for i in range(rl)])
	ppos = 0
	flag = 0
	
	sys.stderr.write("hashing fasta index...\n")
	hfai,flast = hash_fai(args.fasta + ".fai")
	
	if hfai is None or flast is None:
		sys.stderr.write("Error: failed to parse fasta index file. Are you sure it exists?\n")
		return 1

	if args.bam_out:
		# make the bam header
		sys.stderr.write("making bam header...\n")
		bheader = make_bam_header(args.fasta + ".fai")	
		# open output file
		sout = pysam.Samfile(args.stub + ".bam", "wb", header = bheader)
	else:
		if args.paired:
			fout1 = open(args.stub + ".1.fq", "w")
			fout2 = open(args.stub + ".2.fq", "w")
		else:
			fout = open(args.stub + ".fq", "w")
		
	fain = open(args.fasta, "r")
	sys.stderr.write("generating reads...\n")
	
	while(rcount < args.total_reads):
		
		# pick random byte offset within fasta file
		offset = int(floor(flast * rand()))
		
		# see if this offset hits any references
		ll = find_loc(hfai, offset)
		
		if len(ll) > 0:
			# got one, let's do this
			
			# determine location within the sequence
			soffset = actual_ppos(offset-ll[0]+1, int(ll[5]), int(ll[6]))
			
			# set fragment length
			if args.paired:
				fl = 0
				while fl < rl:
					# fl = int(floor(gauss(args.f_mean, args.f_sdev)))
					fl = int(floor(normal(frag_len, frag_sdev)))
			else:
				fl = rl
			
			if soffset + fl < int(ll[4]):
				
				# fragment length fits within the sequence at this offset - go!
				
				# make a read
				# jump to location then read rl*2 bases either forward or backward
				# strip out newlines
				# trim read down to read length
				# mutate bases, etc
				# generate read name
				# print it out
				
				# if we land at 456 with read length 100 and want to read the reverse
				# strand then we need to jump back to 456-100+1 = 357
				
				fain.seek(offset)
				rbuff = fain.read(fl*2).upper()
				# strip out newlines
				rbuff = rbuff.replace("\n", "")
				rbuff = rbuff[0:fl]
				
				# check to make sure that we only have sequence in there
				m = re.search("[^ACTGNU]", rbuff)
				if (not m) and (len(rbuff)==fl):
					
					rev = rand() < 0.5

					query = []
					errors = []
					rmm = []
					
					if args.paired:
						# first read
						query.append(rbuff[0:rl])
						# second read
						query.append(reverse_comp(rbuff[(len(rbuff)-rl):]))
						
						# make errors
						res = generate_base_errors(query[0], args.e_rate)
						query[0] = res[0]
						if len(res[1]) == 0:
							errors.append("NE")
							rmm.append(0)
						else:
							errors.append(res[1])
							rmm.append(res[1].count(">"))

						res = generate_base_errors(query[1], args.e_rate)
						query[1] = res[0]
						if len(res[1]) == 0:
							errors.append("NE")
							rmm.append(0)
						else:
							errors.append(res[1])
							rmm.append(res[1].count(">"))
						
						qname = "{:s}_{:d}_{:d}_{:d}_{:d}_{:s}_{:s}".format(
								ll[2], soffset, soffset-1+fl, int(rev), rcount, errors[0], errors[1])
						
					else:
						
						query.append(rbuff)
						
						# introduce errors
						res = generate_base_errors(query[0], args.e_rate)
						query[0] = res[0]
						
						if len(res[1]) == 0:
							errors.append("NE")
							rmm.append(0)
						else:
							errors.append(res[1])
							rmm.append(res[1].count(">"))
												
						# make read name
						qname = "{:s}_{:d}_{:d}_{:d}_{:d}_{:s}".format(ll[2], soffset, soffset-1+rl, int(rev), rcount, errors[0])


					if args.bam_out:
						a = pysam.AlignedRead()
					
						if args.paired:
							a1 = pysam.AlignedRead()
							a2 = pysam.AlignedRead()

							a1.flag = 0x40+0x2+0x1
							a2.flag = 0x80+0x2+0x1
							
							
							#
							# if read is reversed then the first read becomes the second
							#
							if rev:
								a1.pos = soffset+(fl-rl)-1
								a1.flag |= 0x10								
								a1.seq = reverse_comp(query[1])
								a1.tags = ( ("XM", rmm[1]), ("NM", rmm[1]) )

								a2.pos = soffset-1
								a2.flag |= 0x20
								a2.seq = query[0]
								a2.tags = ( ("XM", rmm[0]), ("NM", rmm[0]) )
								
								
							else:
								a1.seq = query[0]
								a1.pos = soffset-1
								a1.flag |= 0x20
								a1.tags = ( ("XM", rmm[0]), ("NM", rmm[0]) )
								
								a2.seq = reverse_comp(query[1])
								a2.pos = soffset+(fl-rl)-1
								a2.flag |= 0x10
								a2.tags = ( ("XM", rmm[1]), ("NM", rmm[1]) )

							a1.qname = qname
							a1.isize = fl
							a1.rname = ll[-1]
							a1.mrnm = a1.rname
							a1.mapq = 255
							a1.cigar = [(0, rl)]
							a1.qual = quals
							a1.mpos = a2.pos
							
							a2.qname = qname
							a2.isize = a1.isize
							a2.rname = a1.rname
							a2.mrnm = a1.rname
							a2.mapq = a1.mapq
							a2.cigar = [(0, rl)]
							a2.mpos = a1.pos
							a2.qual = quals							
							
							sout.write(a1)
							sout.write(a2)
						
						else:
							a = pysam.AlignedRead()
							
							a.flag = 0
							if rev:
								# reverse compliment
								# rbuff = reverse_comp(rbuff)
								a.flag |= 0x10
						
							a.qname = qname
							a.seq = res[0]
							a.rname = ll[-1]
							a.pos = soffset-1
							a.mapq = 255
							a.cigar = [(0, rl)]
							a.mrnm = -1
							a.mpos = -1
							a.isize = 0
							a.qual = quals
							a.tags = ( ("XM", rmm[0]), ("NM", rmm[0]) )
							
							sout.write(a)
						
					else:
						
						if args.paired:
							
							if rev:
								# have to switch the reads if reversed but no need to reverse compliment
								q0 = query[0]
								query[0] = query[1]
								query[1] = q0
							
							# write read to files
							
							fout1.write("@" + qname + "/1\n")
							fout1.write(query[0] + "\n")
							fout1.write("+\n")
							fout1.write(quals + "\n")

							fout2.write("@" + qname + "/2\n")
							fout2.write(query[1] + "\n")
							fout2.write("+\n")
							fout2.write(quals + "\n")
						
						else:
							
							if rev:
								query[0] = reverse_comp(query[0])
							
							# write read to file
							fout.write("@" + qname + "\n")
							fout.write(query[0] + "\n")
							fout.write("+\n")
							fout.write(quals + "\n")
							
						
					rcount += 1


	fain.close()
	
	if args.bam_out:
		sout.close()
	else:
		if args.paired:
			fout1.close()
			fout2.close()
		else:
			fout.close()

	return 0

#===============================================================================
# hash_fai
# reads the FASTA index into a quick lookup hash where references in the 
# index are hashed by their location in the fasta file.
#===============================================================================
def hash_fai(fname):
	
	fai = {}
	sstart = 0
	slast = 0
	rid_hash = ""
	ll = []
	szl = ""
	ifl = 0
	ref_index = 0
	
	# open fai file
	try:
		fin = open(fname, "r")
	except IOError,e:
		return None,None
	
	for szl in fin:
		ll = szl.strip().split("\t")
		sstart = int(ll[2])
		ifl = seq_len(ll)
		slast = sstart + ifl - 1
		
		rid_hash = hash_loc(sstart)
		
		if rid_hash not in fai:
			fai[rid_hash] = []
		
		fai[rid_hash].append([sstart, slast]+list(ll) + [ref_index])
		ref_index += 1
	
	fin.close()
	
	return fai,slast

def find_loc(fai, loc):
	
	res = []
	
	# hash the location
	loc_hash = hash_loc(loc)
	
	if loc_hash not in fai:
		return []
	
	ll = fai[loc_hash]
	for i in range(len(ll)):
		if loc >= ll[i][0] and loc <= ll[i][1]:
			res = list(ll[i])
			break
	
	return res
		

def hash_loc(n):
	nbin = 16000
	return hashlib.md5(str(floor(n/nbin))).hexdigest()

def seq_len(ll):
	
	# number of lines less any partial line at the end
	n1 = int(ll[1])/int(ll[3])
	# total length in file
	n2 = (n1*int(ll[4])) + (int(ll[1]) - n1*int(ll[3]))
	
	return n2

#===============================================================================
# actual_ppos
# translate offset within a reference that includes newlines to actual offset
# without newlines.
#===============================================================================
def actual_ppos(ppos, n1, n2):
	
	# find number of times the position is divisible by the newline line length
	# lnum = int(floor(ppos/n2))
	
	nppos = ppos - int(floor(ppos/n2))
	if nppos % n1 == 0:
		nppos += 1
	
	return nppos
	
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
	prob = rand(read_length)
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

	return list([read, e_string])

def make_bam_header(fname):
	
	header = {}
	
	try:
		fin = open(fname, "r")
	except IOError, e:
		return 1
	
	header = { 'HD': {'VN': '1.0'}, 'SQ': [] }
	
	for szl in fin:
		ll = szl.strip().split("\t")
		header['SQ'].append({'LN': int(ll[1]), 'SN': ll[0]})
	
	fin.close()
	return header

#===============================================================================
# main entry point, parse args and call main
#===============================================================================

parser = argparse.ArgumentParser(description="Short read (SE or PE) simulator for references.")
parser.add_argument('fasta', type=str, help="FASTA input file")
parser.add_argument('-N', type=int, default=1000000, dest="total_reads", 
				help="Total reads to generate (default: 1000000)")
parser.add_argument('-e', type=float, default=0.01, dest="e_rate", 
				help="Base-call error rate (default: 0.01)")
parser.add_argument('-o', dest="stub", type=str, default="reads", help="Output stub for output files")
parser.add_argument('--bam', dest="bam_out", action="store_const", const=True, default=False,
				help="Output reads as bam alignments (default: off)")
parser.add_argument('-l', dest="read_length", action="store", type=int, default=50, 
				help="Read length in bp (default: 50)")
parser.add_argument('--paired', dest="paired", action="store_const", const=True, default=False,
				help="Generate paired-end reads (default: off)")
parser.add_argument('-fm', dest="f_mean", action="store", default=250, type=int, 
				help="Mean fragment length for paired end data, with --paired (default: 250")
parser.add_argument('-fs', dest="f_sdev", action="store", type=int, default=80, 
				help="Standard deviation of fragment length, with --paired (default: 80)")


args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
	

