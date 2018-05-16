#!/usr/bin/python
#
# Shawn Driscoll
# 20170510
#
# Procedural style SAM parsing. Maybe 25% faster than OOP.
#

import sys, re
import subprocess as sp

#
# parse a SAM record into a dict. if the read is aligned then the boundaries
# of the alignment are calculated as well as the aligned length
def samaln_parse(sz):
	aln = {}

	tmp = sz.strip().split("\t")

	aln['qname'] = tmp[0]
	aln['flag'] = int(tmp[1])
	aln['rname'] = tmp[2]
	aln['pos'] = int(tmp[3])
	aln['mapq'] = int(tmp[4])
	aln['cigar'] = tmp[5]
	aln['rnext'] = tmp[6]
	aln['pnext'] = int(tmp[7])
	aln['tlen'] = int(tmp[8])
	aln['seq'] = tmp[9]
	aln['qual'] = tmp[10]

	aln['read_len'] = len(tmp[9])

	# find boundaries of the alignment as well as the aligned length which 
	# would be different from the read length in the event of soft-clipping
	if not samaln_unaligned(aln):
		samaln_calc_bounds(aln)

	aln['attr'] = {}

	# parse out attributes
	if len(tmp) > 11:
		for i in range(11, len(tmp)):
			tparts = tmp[i].split(":")
			value = tparts[2]
			if tparts[1] == "i":
				value = int(value)
			elif tparts[1] == "f":
				value = float(value)

			aln['attr'][tparts[0]] = value

	return aln

#
# functions to check the status of the alignment flag
def samaln_paired(aln):
	return True if (aln['flag'] & 0x1) != 0 else False

def samaln_properly_paired(aln):
	return True if (aln['flag'] & 0x1) != 0 else False

def samaln_unaligned(aln):
	return True if (aln['flag'] & 0x4) != 0 else False

def samaln_mate_unaligned(aln):
	return True if (aln['flag'] & 0x8) != 0 else False

def samaln_is_reversed(aln):
	return True if (aln['flag'] & 0x10) != 0 else False

def samaln_mate_reversed(aln):
	return True if (aln['flag'] & 0x20) != 0 else False

def samaln_is_first_mate(aln):
	return True if (aln['flag'] & 0x40) != 0 else False

def samaln_is_second_mate(aln):
	return True if (aln['flag'] & 0x80) != 0 else False

def samaln_is_secondary(aln):
	return True if (aln['flag'] & 0x100) != 0 else False

#
# compute the boundaries of an alignment
def samaln_calc_bounds(aln):
	aln['bounds'] = []

	left = right = aln['pos']

	op_len = re.split("[A-Z]", aln['cigar'])
	op_typ = re.split("[0-9]+", aln['cigar'])
	
	# drop the extra bit that comes along with the split
	op_len = map(int, op_len[0:(len(op_len)-1)])
	op_typ = op_typ[1:len(op_typ)]
	nop = len(op_typ)

	# loop through
	aln['aln_len'] = 0
	for i in range(nop):
		if op_typ[i] == "M" or op_typ[i] == "D":
			right += op_len[i]
			if op_typ[i] == "M":
				aln['aln_len'] += op_len[i]

		if op_typ[i] == "N":
			# start a new feature and pass the current on
			aln['bounds'].append([left, right-1])
			left = right + op_len[i]
			right = left

	# append feature
	aln['bounds'].append([left, right-1])

	return 0

if __name__ == "__main__":

	# assume bam file passed at the command line
	if len(sys.argv) < 2:
		sys.stderr.write("useage: {} <input.bam>\n".format(sys.argv[0]))
		sys.exit(1)

	bamfile = sys.argv[1]

	# to read in a bam we need to call on samtools
	cmd = "samtools view {}".format(bamfile)
	p1 = sp.Popen(cmd.split(), stdout=sp.PIPE)
	fin = p1.stdout

	# read lines, count the number that are aligned and count number in each MAPQ bin
	mapq = {}
	num_parsed = num_aligned = 0
	for szl in fin:
		aln = samaln_parse(szl)
		num_parsed += 1
		if not samaln_unaligned(aln):
			num_aligned += 1
			if aln['mapq'] not in mapq:
				mapq[aln['mapq']] = 0
			mapq[aln['mapq']] += 1

		# ... do other stuff 

	print "parsed {} lines; found {} aligned reads".format(num_parsed, num_aligned)

	print "MAPQ histogram"
	for k in mapq.keys():
		print "{}\t{}".format(k, mapq[k])

	# no need to close connection

	sys.exit(0)
	


