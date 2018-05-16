#!/usr/bin/env python
#
# tsam-to-gsam.py
#
# Shawn Driscoll
# 20121206
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script translates transcriptome alignments to genome alignments in SAM 
# format. A GTF file is required that describes the transcriptome elements
# in genome coordinates. The reference id's in the SAM data must match
# the transcript names in the GTF file.
#

import sys, argparse, re
import subprocess as sp
import pysam

#
# main script
#

def main(args):
	#
	# variables
	#

	output_buffer = []
	paired_ppos_buffer =[]
	paired_ref_buffer = []
	paired = False
	input_buffer = []
	min_anchor = 2

	#
	# check arguments. if the second argument is missing it is assumed that the same
	# input is passed via stdin and the first argument is the GTF file
	#
	use_stdin = False
	if args.sam == "-":
		use_stdin = True
	else:
		sam = args.sam
	
	# grab name of gtf
	gtf = args.gtf
		
	# parse the gtf
	dgtf = parse_gtf(gtf)
	
	if use_stdin:
		fin = sys.stdin
		if args.sam_in:
			fin = pysam.Samfile("-", "r")
	else:
		fin = open(sam, "r")
	
	aln_rev = False
	
	for szl in fin:
		szl = szl.strip()

		# skip header lines
		if len(szl) > 0 and szl[0] == "@":
			pass
		elif len(szl) > 0:
			ll = szl.split("\t")
			flag = int(ll[1])

			# paired?
			paired = flag & 0x1 != 0

			input_buffer.append(list(ll))

			# keep input buffer locked at two alignments
			if len(input_buffer) > 2:
				input_buffer.pop(0)

			# 
			# gather information from the alignment
			#	
			tid = ll[2]
			ppos_ref = int(ll[3])
			types,lengths = parse_cigar(ll[5])
			flag = int(ll[1])
			aln_rev = False
			
			
			if not paired and flag & 0x4:
				# unaligned single-end read, print it back out
				sys.stdout.write("\t".join(map(str, ll)) + "\n")

			else:

				if flag & 0x4:
					# if paired and the current read is unaligned the code will make it
					# into this condition.
					output_buffer.append(list(ll))

				else:

					# is the alignment reversed?
					aln_rev = flag & 0x10

					# kill off the MD tag since i'll likely not convert it properly
					i = 11
					while ll[i][0:2] != "MD" and i < len(ll):
						i += 1
					
					if i < len(ll):
						if ll[i][0:2] == "MD":
							# reverse the MD tag
							#parts = ll[i].split(":")
	
							#ll[i] = "MD:Z:" + reverse_md(parts[2])
							ll.pop(i)
#							if i+1 < len(ll):
#								j = i+1
#								while j < len(ll):
#									ll[i] = ll[j]
#									i += 1
#									j += 1
	
#							temp = list(ll[0:(len(ll)-1)])
#							ll = list(temp)

					
					# get total alignment length
					aln_length = 0
					for i in range(len(types)):
						if types[i] == "M" or types[i] == "D":
							aln_length += lengths[i]
					
					# fetch the feature from the hashed gtf
					if tid not in dgtf:
						sys.stderr.write("Warning: unable to locate transcript id the read was aligned to ({:s})\n".format(tid))
						# set unaligned
						ll[1] = set_bit(flag, 0x4)
						if paired:
							output_buffer.append(list(ll))
						else:
							sys.stdout.write("\t".join(map(str, ll)) + "\n")

					else:

						# gather feature information from GTF hash
						strand = dgtf[tid]['strand']
						glist = dgtf[tid]['g_features']
						tlist = dgtf[tid]['t_features']
						num_exons = len(glist)
						direction = strand == "+"
											
						# fix the flag and the read. if the isoform is on the - strand then
						# everything needs to be reversed
						if not direction:
							if aln_rev:
								# ll[1] = str(int(ll[1]) - 16)
								ll[1] = str(clear_bit(flag, 0x10))
							else:
								# ll[1] = str(int(ll[1]) + 16)
								ll[1] = str(set_bit(flag, 0x10))

							# reverse compliment the read bases
							ll[9] = rev_comp(ll[9])
							# reverse the base qual string
							ll[10] = ll[10][::-1]
							# reverse all of the CIGAR data
							types.reverse()
							lengths.reverse()							

						#if
						
						# first figure out which exon this alignment starts in in terms of 
						# transcriptome
						found = False
						i = 0
						while i < num_exons and not found:
							if ppos_ref >= tlist[i][0] and ppos_ref <= tlist[i][1]:
								found = True
								
							i += 1

						#while
						
						i -= 1
						start_exon = tlist[i]
						diff = ppos_ref - start_exon[0]
						#print start_exon
						#print diff
						
						# now based on the strand of the feature we can translate the start position				
						if direction:
							g_index = i
							genome_start = glist[g_index][1] + diff
						else:
							# if feature strand is negative then the genome start position is really the
							# transcriptome END position of the alignment which means the transcriptome
							# exon may not be the same genome exon that the alignment starts in, in the case
							# of an alignment that crosses a splice junction
							
							# translate transcriptome exon to genome exon index
							g_index = num_exons - 1 - i
							genome_start = glist[g_index][2] - diff - aln_length + 1
							# check if this means it starts in a different exon in genome
							# coordinates.
							if genome_start < glist[g_index][1]:
								# in some cases a read might span several exons - loop until
								# the start position is contained
								while genome_start < glist[g_index][1] and g_index >= 0:								
									diff = glist[g_index][1] - genome_start 
									g_index -= 1
									genome_start = glist[g_index][2] - diff + 1
						#if
								
						#print g_index
						#print glist[g_index]
						#print direction

						# loop through the CIGAR fields and build new alignment CIGAR
						i = 0
						n = len(types)
						ppos_last = genome_start
						exon_last = g_index
						new_cigar = ""
						drop_alignment = False
						direction = strand == "+"
						
	#					print glist
	#					print tlist
	#					print strand, genome_start, ppos_ref
						
						while i < n and not drop_alignment:
							if types[i] == "M":
								
								# since I realized if the isoform is on the - strand that the alignment
								# start position for the genome is actually the end of the alignment in the
								# transcriptome then by flipping the CIGAR we can work through the alignment
								# in the same way for both + and - strand features

								# find end coordinate
								temp_end = ppos_last + lengths[i] - 1
	
								# did this just jump an exon?
								if temp_end > glist[exon_last][2]:
									# in some cases a single chunk may leap across a couple junctions.
									# loop until the remaining alignment bases are contained
									rem = lengths[i]
									while temp_end > glist[exon_last][2]:
										# compute number of bases that fits into the current exon
										diff = glist[exon_last][2] - ppos_last + 1
										new_cigar += str(diff) + "M"
										# how many bases are left in the alignment?
										rem = rem - diff
										# insert length of intron
										diff = glist[exon_last+1][1] - glist[exon_last][2] - 1
										new_cigar += str(diff) + "N"
										
										# increment exon index
										exon_last += 1
										ppos_last = glist[exon_last][1]
										temp_end = ppos_last + rem - 1
									#while
									
									# temp_end doesn't extend past the current exon so we can report the match
									# insert the remainder of the matched bases
									new_cigar += str(rem) + "M"
									# update position to the next possible base
									ppos_last = glist[exon_last][1] + rem
									
										
								else:
									# alignment is totally contained in this exon
									new_cigar += str(lengths[i]) + "M"
									ppos_last += lengths[i]

								#if
							elif types[i] == "D":
								# for deletion type we need to check if the deletion itself crosses a junction.
								# if it does then it's probably best to drop this alignment
								
								if True:
								
									temp_end = ppos_last + lengths[i] - 1
			
									# did this just jump an exon?
									if temp_end > glist[exon_last][2]:
										drop_alignment = True
									else:
										# deletion doesn't cross a junction so it is just passed on
										# to the new cigar string and the genomic position is updated
										new_cigar += str(lengths[i]) + "D"
										ppos_last += lengths[i]
									#if
								else:
									
									temp_end = ppos_last - (lengths[i] - 1)
			
									# did this just jump an exon?
									if temp_end < glist[exon_last][1]:
										drop_alignment = True
									else:
										# deletion doesn't cross a junction so it is just passed on
										# to the new cigar string and the genomic position is updated
										new_cigar += str(lengths[i]) + "D"
										ppos_last -= lengths[i]
									#if
								#if
							elif types[i] == "I":
								# insertion doesn't change anything in the genomic alignment so we can just
								# put this back into the cigar string
								new_cigar += str(lengths[i]) + "I"
							elif types[i] == "S":
								# soft clipping is just passed back to the cigar and it doesn't 
								# alter the genome position since they are bases ignored at the end or the
								# start of the alignment
								new_cigar += str(lengths[i]) + "S"
							else:
								sys.stderr.write("Warning: unexpected CIGAR value {:s}. Skipping alignment.".format(types[i]))
								drop_alignment = True
							#if
							
							i += 1

						#while
						
						if drop_alignment:

							# set as unaligned - we'll fix up the rest of the flags down the road
							ltmp = set_unaligned(ll)
							#ll[1] = str(set_bit(int(ll[1]), 0x4))
							#ll[2] = "*"
							#ll[5] = "*"

							if paired:
								output_buffer.append(list(ltmp))
							else:
								sys.stdout.write("\t".join(map(str, ltmp)) + "\n")
							#}
							
						else:
							#
							# finalize the new alignment
							#
							
							# the new cigar string should be checked for anchors on one side or the other of a 
							# junction, 'N' region, less than some number of bases. if that's true then
							# instead of dropping the hit the shorter of the two will be converted to a soft-clip
							# event

							new_types, new_lengths = parse_cigar(new_cigar)
							if new_types[0] == "M" and new_lengths[0] < min_anchor and new_types[1] == "N":
								# turn the M into a soft clip and eliminate the N region
								new_types[0] = "S"
								new_types[1] = "*"
								# it's also necessary to update the alignment start position to exclude the
								# junction hit
								genome_start += (int(new_lengths[1]) + int(new_lengths[0]))
							elif new_types[-1] == "M" and new_lengths[-1] < min_anchor and new_types[-2] == "N":
								# turn M into a soft clip and eliminate the N
								new_types[-1] = "S"
								new_types[-2] = "*"
							#}
							# rebuild string from lists
							ltmp = []
							for i in range(len(new_types)):
								if new_types[i] != "*":
									ltmp.append(new_lengths[i])
									ltmp.append(new_types[i])
								#}
							#}

							new_cigar = "".join(map(str, ltmp))

							# change transcript id to genome reference
							ll[2] = glist[0][0]	
							# update alignment position							
							ll[3] = str(genome_start)
							# update CIGAR
							ll[5] = new_cigar
							# update MAPQ - this is obviously meaningless
							ll[4] = str(10)
							

							if paired:
								# buffer the output of paired alignments to provide a chance to correct
								# RNEXT, PNEXT and insertion size fields for pairs
								output_buffer.append(list(ll))
							else:
								sys.stdout.write("\t".join(map(str, ll)) + "\n")
							#}
						#}

					#} end hashed feature exists in GTF

				#} end feature is aligned
				
				# address the output buffer if the data is paired
				if paired:

					if len(output_buffer) > 1:
						result = are_paired(input_buffer[0], input_buffer[1])

						flag1 = int(output_buffer[0][1])
						flag2 = int(output_buffer[1][1])

						if result >= 0:

							# double check the current flag status because the pair's relationship may
							# be different
							if flag1 & 0x4 and flag2 & 0x4:
								# now this is an unaligned pair
								result = 3
							elif flag1 & 0x4 or flag2 & 0x4:
								# now this is a mixed alignment
								result = 2
							#}

							if result == 3:
								# this is an unaligned pair
								output_buffer[0] = list(set_unaligned(output_buffer[0]))
								output_buffer[1] = list(set_unaligned(output_buffer[1]))

							elif result == 2:

								ltmp = set_mixed_pair(output_buffer)
								output_buffer = list(ltmp)

							else:

								# fix the pair info fields
								if output_buffer[0][2] == output_buffer[1][2]:
									output_buffer[0][6] = "="
									output_buffer[1][6] = "="
								else:
									output_buffer[0][6] = output_buffer[1][2]
									output_buffer[1][6] = output_buffer[0][2]
								#}

								output_buffer[0][7] = output_buffer[1][3]
								output_buffer[1][7] = output_buffer[0][3]

								# fix insert size orientation
								isize = abs(int(output_buffer[0][8]))

								if int(output_buffer[0][1]) & 0x10:
									output_buffer[0][8] = str(-1 * isize)
								else:
									output_buffer[0][8] = str(isize)
								#}

								if int(output_buffer[1][1]) & 0x10:
									output_buffer[1][8] = str(-1 * isize)
								else:
									output_buffer[1][8] = str(isize)
								#}

								# confirm that the 0x20 flag is set correctly
								if flag1 & 0x10 and not flag2 & 0x10:
									flag2 = set_bit(flag2, 0x20)
									flag1 = clear_bit(flag1, 0x20)
								elif not flag1 & 0x10 and flag2 & 0x10:
									flag1 = set_bit(flag1, 0x20)
									flag2 = clear_bit(flag2, 0x20)

								output_buffer[0][1] = str(flag1)
								output_buffer[1][1] = str(flag2)

							#}

							sys.stdout.write("\t".join(map(str, output_buffer[0])) + "\n")
							sys.stdout.write("\t".join(map(str, output_buffer[1])) + "\n")
							output_buffer = []
							input_buffer = []
							paired_ppos_buffer = []
							paried_ref_buffer =[]

						else:
							# result < 0, not a pair
							sys.stderr.write("Warning: processing orphaned mate: {:s}\n".format(output_buffer[0][0]))

							sys.stdout.write("\t".join(map(str, output_buffer[0])) + "\n")
							output_buffer.pop(0)
							input_buffer.pop(0)
							paired_ppos_buffer.pop(0)
							paired_ref_buffer.pop(0)

						#}

					#} # buffer not full

				#}	not paired		

			#} # finished with dealing with aligned read or unaligned mate of a pair

		#} # finished with 'not sam header' section

	#} # end of main while loop
	
	if len(output_buffer) > 0:
		sys.stderr.write("Warning: {:d} items remaining in buffer. printing..\n".format(len(output_buffer)))
		for i in range(len(output_buffer)):
			sys.stdout.write("\t".join(map(str, output_buffer[i])) + "\n")
			
	if not use_stdin:
		fin.close()

	# end main

def reverse_md(sz):

	toggle = 0
	i = 0
	n = len(sz)
	parts = []
	while i < n:
		if toggle == 0:
			m = re.search("([0-9]+)", sz[i:])
			res = m.group(1)
			parts.append(res)
			i += len(res)
		elif toggle == 1:
			m = re.search("([A-Z]|\^[A-Z]+)", sz[i:])
			res = m.group(1)

			if res[0] == "^":
				temp = rev_comp(res[1:])
				res = "^" + temp
			else:
				res = rev_comp(res)

			parts.append(res)
			i += len(res)

		toggle = toggle * -1 + 1

	parts.reverse()

	return "".join(parts)

def parse_cigar(sz):
	#
	# parse CIGAR notation into two lists. one list of CIGAR types and another
	# of the lengths associated with the type
	#
	split_1 = re.split("[A-Z]",sz)
	split_2 = re.split("[0-9]+",sz)

	split_1 = split_1[0:(len(split_1)-1)]
	split_2 = split_2[1:]

	return (split_2,map(int,split_1))

def rev_comp(read):
	#
	# reverse compliment a string of ACTGN's
	#
	ttable = {"A":"T", "C":"G", "T":"A", "G":"C", "N":"N"}
	n = len(read)
	rev_read = ""
	while n > 0:
		n -= 1
		rev_read += ttable[read[n]]

	return rev_read

#
# parse_gtf
# Parses GTF annotation (exons only) to build a translation table for alignments directly
# to the transcripts to alignments in genome coordinates
def parse_gtf(gtf):
	# sort the exon list within the transcript
	# dgtf[tid]['features'].sort(key=lambda x: int(x[3]))
	
	type = "exon"
	
	# open file
	try:
		fin = open(gtf, "r")
	except IOError, e:
		sys.stderr.write("Error opening input GTF file {:s}\n".format(gtf))
		return None
	
	# onward...
	
	dgtf = {}

	for szl in fin:

		ll = szl.strip().split("\t")
		if ll[2] == type:
			attr = parse_gtf_attr(ll[8])
	
			tid = attr['transcript_id']
	
			# if transcript id isn't present in the locus, add it
			if tid not in dgtf:
				dgtf[tid] = {}
				dgtf[tid]['g_features'] = []
				dgtf[tid]['strand'] = ll[6]
	
			# append feature row list to transcript within locus
			dgtf[tid]['g_features'].append([ll[0], int(ll[3]), int(ll[4]) ])
			
	fin.close()

	# make sure that for each transcript list the exon features are sorted by start position. 
	# once they are sorted the transcriptome coordinates can be generated
	
	for tid in dgtf.keys():
		# tid is current transcript id
		
		# sort the exon list within the transcript
		dgtf[tid]['g_features'].sort(key=lambda x: int(x[1]))
		
		#
		# produce transcriptome coordinates
		#
		
		# get exon list
		elist = dgtf[tid]['g_features']
		# initalize a new list for the transcriptome features
		dgtf[tid]['t_features'] = []
		
		temp = []
		n = len(elist)
		
		# loop from start to end for + strand features and from end to start 
		# for - strand features
		
		if dgtf[tid]['strand'] == "+":
			i = 0
			e_start = 1
			while i < n:
				exon_length = elist[i][2] - elist[i][1] + 1
				temp.append([e_start, e_start + exon_length - 1])
				e_start += (exon_length)
				i += 1
				
			dgtf[tid]['t_features'] = list(temp)
		else:
			i = n
			e_start = 1
			while i > 0:
				i -= 1
				exon_length = elist[i][2] - elist[i][1] + 1
				temp.append([e_start, e_start + exon_length - 1])
				e_start += (exon_length)

			dgtf[tid]['t_features'] = list(temp)
				
	return(dgtf)

# 
# parse_gtf_attr
def parse_gtf_attr(field):

	attrs = {}

	# split into fields
	l1 = field.split(";")

	# split each field into key and value
	l2 = []
	for i in range(len(l1)):
		l2.append(l1[i].split("\""))


	for i in range(len(l2)):
		if len(l2[i]) > 1:
			attrs[l2[i][0].strip()] = l2[i][1].strip()

	return(attrs)		
	
def are_paired(r1, r2):

	flag1 = int(r1[1])
	flag2 = int(r2[1])

	unaligned = flag1 & 0x4 and flag2 & 0x4
	possible_mixed = (flag1 & 0x8 and flag2 & 0x4) or (flag1 & 0x4 and flag2 & 0x8)
	fns = (flag1 & 0x40 and flag2 & 0x80) or (flag1 & 0x8 and flag2 & 0x40)
	possible_mates = r1[0] == r2[0] and fns
	possible_discordant = r1[6] != "=" or r2[6] != "="

	if r1[0] != r2[0]:
		# diffferent read names
		return -1
	#}
	
	if flag1 & 0x4 and flag2 & 0x4:
		# unaligned pair
		return 3
	#}
	
	if r1[6] != "=" and r2[6] != "=":
		# possible discordant
		if r1[6] == r2[2] and r2[6] == r1[2]:
			return 1
		else:
			return -1
		#}
	#}
	
	# check for singleton - this is a little confused now because if a read
	# is discarded it is simply set to unaligned. in the paired data then that
	# read's mate may still be aligned and we'll have to just trust that since
	# the read names are equal and it isn't a discordant alignment that these
	# are a pair
	
	# if (flag1 & 0x4 and flag2 & 0x8) or (flag1 & 0x8 and flag2 & 0x4):
	if (flag1 & 0x4) or (flag2 & 0x4):
		# possible mixed
		
		# the only thing I can think of to check is if the matched read names are't
		# both primary or both secondary alignments. this entire code file is trusting that
		# the data is coming in sorted with pairs right next to each other so this test is 
		# kinda pointless...
		if ((flag1 & 0x100) and not (flag2 & 0x100)) or ((flag2 & 0x100) and not (flag1 & 0x100)):
			return -1
		
		return 2
	
#		if int(r1[3]) == int(r2[3]):
#			return 2
#		else:
#			return -1
		#}
	#}
	
	# one last thing to check
	if int(r1[7]) != int(r2[3]) or int(r2[7]) != int(r1[3]):
		return -1
	#}

	if False:
		if not possible_mates:
			return -1
	
		if unaligned:
			return 3
	
		if possible_mixed:
			# confirm
			if r1[3] == r2[3] and r1[7] == r2[7] and r1[2] == r2[2]:
				return 2
			else:
				return -1
	
		if possible_discordant:
			# confirm
			if r1[6] == r2[2] and r2[6] == r1[2] and r1[7] == r2[3] and r1[3] == r2[7]:
				return 1
	
			else:
				return -1
	
		# confirm proper alignment
		if r1[6] != "=" or r2[6] != "=":
			return -1
	
		if r1[7] != r2[3] or r1[3] != r2[7]:
			return -1

	return 0
#}

def set_bit(v, bit):
	return v | bit
#}

def clear_bit(v, bit):
	return v & ~bit
#}

def set_unaligned(aln):
	# aln is a list of a parsed sam alignment

	flag = int(aln[1])
	new_flag = 0x1 + 0x4
	aln_new = []
	ext = ""

	if flag & 0x1:
		if flag & 0x40:
			new_flag += 0x40
		else:
			new_flag += 0x80
		#}
		ext = "YT:Z:UP"
	else:
		new_flag = 0x4
		ext = "YT:Z:UU"
	#}

	aln_new = [aln[0], new_flag, "*", 0, 0, "*", "*", 0, 0, aln[9], aln[10], ext]

	return aln_new
#}


def set_mixed_pair(lbuffer):

	r1 = list(lbuffer[0])
	r2 = list(lbuffer[1])
	flag1 = int(r1[1])
	flag2 = int(r2[1])
	
	aligned = 1
	unaligned = 0
	if flag2 & 0x4:
		aligned = 0
		unaligned = 1
		flag1 = set_bit(flag1, 0x8)
	else:
		flag2 = set_bit(flag2 , 0x8)
	#}
	
	aligned_read = lbuffer[aligned]
	
	# zero insert size
	r1[8] = 0
	r2[8] = 0
	# fix RNEXT and PNEXT fields
	r1[6] = "="
	r2[6] = "="
	r1[7] = aligned_read[3]
	r2[7] = aligned_read[3]

	# fix CIGAR of unaligned mate
	if aligned == 0:
		r2[5] = "*"
		r2[2] = r1[2]
		r2[3] = r1[3]
	else:
		r1[5] = "*"
		r1[2] = r2[2]
		r1[3] = r2[3]
	#}
	
	# make sure the flags make sense
	
	if aligned == 0:
		if flag1 & 0x10:
			flag2 = set_bit(flag2, 0x20)
		else:
			flag2 = clear_bit(flag2, 0x20)
		#}
	else:
		if flag2 & 0x10:
			flag1 = set_bit(flag1, 0x20)
		else:
			flag1 = clear_bit(flag1, 0x20)
		#}
	#}
	
	r1[1] = str(flag1)
	r2[1] = str(flag2)
	
	return [r1, r2]
	
#}

#==============================================================================
# main entry point
#==============================================================================

#
# parse arguments
#

parser = argparse.ArgumentParser(description="Translate transcripome aligned alignments to genome coordiantes.")
parser.add_argument("bam", type=str, help="Alignments in BAM format or - to read from STDIN")
parser.add_argument("gtf", type=str, help="Gene annotation in GTF format")
parser.add_argument("-S", dest="sam_in", action="store_const", const=True, default=False, 
				help="Input is SAM and not BAM")


args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))








