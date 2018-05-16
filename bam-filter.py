#!/usr/bin/env python
#==============================================================================
# bam-filter.py
#
# Shawn Driscoll
# 20130409
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Attempt to use pysam to parse a BAM file and filter based on user options.
# Default (without options) produces some simple stats.
#==============================================================================

import sys, re, argparse
import pysam

#
# main()
def main(args):
	#
	# variables
	#
	sf = None
	fin = None
	i = 0
	aln = None
	paired = False
	aln_buffer = []
	processed = []
	result = 0
	cnt = 0

	#
	# open bam file
	#
	if args.S:
		if args.bam == "-" or args.bam == "stdin":
			sf = pysam.Samfile("-", "r")
		else:
			sf = pysam.Samfile(args.bam, "r")
	else:
		sf = pysam.Samfile(args.bam, "rb")
	#}

	# I guess output will just be in the same format as the input
	if len(args.o) > 0:
		if args.s:
			samout = pysam.Samfile(args.o, "w", template=sf)
		else:
			samout = pysam.Samfile(args.o, "wb", template=sf)
	else:
		if args.s:
			samout = pysam.Samfile("-", "w", template=sf)
		else:
			samout = pysam.Samfile("-", "wb", template=sf)
	#}
	
	for aln in sf:

		cnt += 1

		if aln.flag & 0x1:
			paired = True

		aln_buffer.append(aln)

		if len(aln_buffer) == 2:
			if paired:
				# paired data

				result = are_paired(aln_buffer[0], aln_buffer[1])

				if result >= 0:

					processed = process_paired_alignment(args, aln_buffer, result)

					aln_buffer = []

					for i in range(2):
						samout.write(processed[i])
					#}
				else:
					# not paired, process orphan read
					sys.stderr.write("Warning: processing orphaned mate: {:s}\n".format(aln_buffer[0].qname))
					processed = process_orphan(args, aln_buffer[0])
					aln_buffer.pop(0)
					samout.write(processed[0])
				#}

			else:
				# process single end reads
				for i in range(2):
					processed = process_alignment(args, aln_buffer[i])
					samout.write(processed[0])
				#}
				aln_buffer = []

			#}

		#}

		if cnt % 1000000 == 0:
			sys.stderr.write("> processed {:d} alignments\n".format(cnt))
		#}

	#}

	sys.stderr.write("> processed {:d} alignments\n".format(cnt))

	sf.close()
	samout.close()

#}


#==============================================================================
# defs
#==============================================================================

def process_paired_alignment(args, lbuffer, result):

	aln1 = lbuffer[0]
	aln2 = lbuffer[1]
	processed = []

	aln1_new = aln1
	aln2_new = aln2
	mm_limit = args.limit_mismatches

	if mm_limit >=0 and mm_limit < 1:
		mm_limit = round(mm_limit * lbuffer[0].rlen)
	#}

	if result < 3:

		if (args.rem_discordant and result == 1) or (args.rem_singletons and result == 2):
			# update pair to unaligned
			aln1_new = set_unaligned(aln1)
			aln2_new = set_unaligned(aln2)
		elif mm_limit >= 0:
			mm1 = get_mismatches(aln1)
			mm2 = get_mismatches(aln2)

			if mm1 > mm_limit or mm2 > mm_limit:
				if args.rem_singletons:
					# don't allow singletons so drop both mates
					aln1_new = set_unaligned(aln1)
					aln2_new = set_unaligned(aln2)
				else:
					# singletons are allowed so in this case the alignments have to be
					# modified to be a mixed alignment
					if mm1 > mm2:
						# drop aln1

						# aln2 becomes the main alignment
						aln2.flag = set_bit(aln2.flag, 0x8)
						aln2.pnext = aln2.pos
						aln2.rnext = aln2.tid
						aln2.tlen = 0

						aln1 = set_unaligned(aln1)
						aln1.pnext = aln2.pos
						aln1.rnext = aln2.tid
						aln1.pos = aln2.pos

					else:
						# aln2 becomes the main alignment
						aln1.flag = set_bit(aln1.flag, 0x8)
						aln1.pnext = aln1.pos
						aln1.rnext = aln1.tid
						aln1.tlen = 0

						aln2 = set_unaligned(aln2)
						aln2.pnext = aln1.pos
						aln2.rnext = aln1.tid
						aln2.pos = aln1.pos
					#}
				#}

			#}
		#}

	#}

	processed.append(aln1_new)
	processed.append(aln2_new)

	return processed
#}

def process_orphan(args, aln):

	aln_new = aln

	mm_limit = args.limit_mismatches

	if mm_limit >=0 and mm_limit < 1:
		mm_limit = round(mm_limit * aln.rlen)
	#}

	if args.rem_singletons:
		aln_new = set_unaligned(aln)
	elif args.rem_secondary and aln.flag & 0x100:
		aln_new = set_unaligned(aln)
	elif mm_limit >= 0:
		mm = get_mismatches(aln)

		if mm > mm_limit:
			aln_new = set_unaligned(aln)
		#}
	#}

	return [aln_new]
#}

def process_alignment(args, aln):

	aln_new = aln

	mm_limit = args.limit_mismatches

	if mm_limit >=0 and mm_limit < 1:
		mm_limit = round(mm_limit * aln.rlen)
	#}

	if args.rem_secondary and aln.flag & 0x100:
		aln_new = set_unaligned(aln)
	elif mm_limit >= 0:
		mm = get_mismatches(aln)

		if mm > mm_limit:
			aln_new = set_unaligned(aln)
		#}
	#}

	return [aln_new]
#}

def set_unaligned(aln):

	if aln.flag & 0x1:
		# paired
		if aln.flag & 0x40:
			aln.flag = 0x1 + 0x4 + 0x40
		else:
			aln.flag = 0x1 + 0x4 + 0x80
		#}
	else:
		aln.flag = 0x4
	#}

	aln.tid = -1
	aln.pos = 0
	aln.pnext = 0
	aln.rnext = -1
	aln.mapq = 0
	aln.cigar = []
	aln.tlen = 0
	aln.tags = ()

	return aln
#}

def get_mismatches(aln):

	# [i for i, v, in enumerate(aln.tags) if v[0] == "FF"]
	indel = 0
	xm = 0
	nm = 0

	total_mm = 0

	# check for xm and nm tags
	ci = get_tag_index(aln.tags, "XM")
	if ci >= 0:
		xm = aln.tags[ci][1]

	ci = get_tag_index(aln.tags, "NM")
	if ci >= 0:
		nm = aln.tags[ci][1]

	if aln.cigar:
		# check cigar for indels
		for i in range(len(aln.cigar)):
			if aln.cigar[i][0] == 1 or aln.cigar[i][0] == 2:
				indel += aln.cigar[i][1]
			#}
		#}
	#}

	total_mm = xm + indel
	if total_mm > nm:
		return total_mm
	#}

	return nm

#}

def get_tag_index(tags, tname):

	result = [i for i, v in enumerate(tags) if v[0] == tname]
	if len(result) == 1:
		return result[0]

	return -1

#}


#}

def file_exists(fname):
	try:
		fin = open(fname)
	except IOError,e:
		return 1

	fin.close()
	return 0
#}

def are_paired(r1, r2):

	flag1 = r1.flag
	flag2 = r2.flag

	unaligned = flag1 & 0x4 != 0 and flag2 & 0x4 != 0
	possible_mixed = (flag1 & 0x8 != 0 and flag2 & 0x4 != 0) or (flag1 & 0x4  != 0 and flag2 & 0x8  != 0)
	fns = (flag1 & 0x40 != 0 and flag2 & 0x80 != 0) or (flag1 & 0x80 != 0 and flag2 & 0x40 != 0)
	possible_mates = r1.qname == r2.qname and fns
	possible_discordant = r1.rnext != r2.rnext

	if not possible_mates:
		return -1

	if unaligned:
		return 3

	if possible_mixed:
		# confirm
		if r1.pos == r2.pos and r1.pnext == r2.pnext:
			return 2
		else:
			return -1

	if possible_discordant:
		# confirm
		if r1.rnext == r2.tid and r2.rnext == r1.tid and r1.pnext == r2.pos and r1.pos == r2.pnext:
			return 1

		else:
			return -1

	# confirm proper alignment
#	if r1.rnext != "=" or r2.rnext != "=":
#		return -1

	if r1.pnext != r2.pos or r1.pos != r2.pnext:
		return -1

	return 0

def set_bit(v, bit):
	return v | bit

def clear_bit(v, bit):
	return v & ~bit


#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Filter alignments in BAM format.")
parser.add_argument("bam", type=str, help="Alignments in BAM format.")
parser.add_argument("-o", dest="o", action="store", type=str, default="",
	help="Output file (default: stdout)")
parser.add_argument("-S", dest="S", action="store_const", const=True, default=False,
	help="Input is SAM format")
parser.add_argument("-s", dest="s", action="store_const", const=True, default=False,
	help="Output SAM instead of BAM (default: off)")
parser.add_argument("--rem-secondary", dest="rem_secondary", action="store_const", const=True,
	default=False, help="Set secondary alignments to unaligned")
parser.add_argument("--rem-singletons", dest="rem_singletons", action="store_const", const=True, 
	default=False, help="Set singletons as unaligned.")
parser.add_argument("--rem-discordant", dest="rem_discordant", action="store_const", const=True, 
	default=False, help="Set discordant pairs as unaligned.")
parser.add_argument("--limit-mismatches", dest="limit_mismatches", type=float, default=-1, action="store", 
	help="Limit number of mismatches per alignment to this value (or ratio of read length")
parser.add_argument("--maximum-insert-size", dest="maximum_insert_size", type=int, default=-1, action="store",
	help="Drop paired alignments with insert size longer than this value")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))

