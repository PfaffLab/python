#!/usr/bin/env python
#==============================================================================
# bam-to-junctions.py
#
# Shawn Driscoll
# 20130626
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script parses BAM or SAM alignments to find junctions (alignments with 
# 'N' cigar operations.  The junctions are reported in BED and a basic 
# info format.  The intron motifs are extracted from the genome so that
# canonical and non-canonical junctions are identified. the motifs are included
# in the .junc file. The motifs are also used to determine the strand of the 
# junction. 
#==============================================================================

import sys, argparse, re, pysam

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	junction_lib = {}
	ppos = 0
	aln_idx = 0
	cu_match_length = 0
	kid = ""
	kid_last = ""
	junctions_in_aln = 0
	dfai = {}
	
	jinfo = args.o + ".junc"
	jbed = args.o + ".bed"
	
	# GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT GTAT
	# jmotifs = ["GTAG", "CTAC", "GCAG", "CTGC", "ATAC", "GTAT"]
	# jmotifs = ["CTAC", "CTGC", "GTAT"]
	jmotifs = ["GTAG", "GCAG", "ATAC"]
	jmotifs_rev = ["CTAC", "CTGC", "GTAT"]
	
	# check input files
	if not file_exists(args.alignments):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.alignments)
		return 1

	if not file_exists(args.fasta):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.fasta)
		return 1
	
	if not file_exists(args.fasta + ".fai"):
		sys.stderr.write("[main] Error: Cannot find .fai index to your fasta file. Please build one with samtools faidx (%s)\n")
		return 1

	# parse the FAI file into a hash
	fin = open(args.fasta + ".fai", "r")
	for szl in fin:
		ll = szl.strip().split("\t")		
		dfai[ll[0]] = map(int, ll[1:])
	
	fin.close()

	# open the alignments as SAM or BAM
	if args.S:
		samin = pysam.Samfile(args.alignments, "r")
	else:
		samin = pysam.Samfile(args.alignments, "rb")
	
	# open the fasta
	fin = open(args.fasta, "r")
	
	# loop through alignments and collect all of the junctions from alignments.
	for aln in samin:
		
		if aln.flag & 0x4:
			# -- skip unaligned --
			pass
		else:
			# -- process aligned reads
			
			aln_idx += 1
			ppos = aln.pos
			
			# progress/count variables for within this alignment
			junctions_in_aln = 0
			cu_match_length = 0
			
			for i in range(len(aln.cigar)):
				
				if aln.cigar[i][0] == 0:
					# track number of aligned bases leading up to the junction
					cu_match_length += aln.cigar[i][1]
					
				if aln.cigar[i][0] == 3:
					# -- found junction
					
					# increment junction count wihin this read
					junctions_in_aln += 1
					
					# get target id and figure out the coordinates of the two ends of the junction
					tid = samin.header['SQ'][aln.tid]['SN']
					el_junction =  [tid, ppos+1, ppos+aln.cigar[i][1]]
					kid = "".join(map(str, el_junction))
					
					if kid not in junction_lib:
						junction_lib[kid] = list(el_junction)
						junction_lib[kid] += [0, 0, 0, ".", ".", 0]
						
						# fetch motif
						jlen = el_junction[2]-el_junction[1]+1
						if jlen >= 4:
							lseq = fetch_seq(fin, dfai, el_junction[0], el_junction[1], el_junction[2])
							lseq = list(lseq.upper())
							jmotif = "".join([lseq[0], lseq[1], lseq[-2], lseq[-1]])

							if jmotif in set(jmotifs):
								junction_lib[kid][6] = "+"
								junction_lib[kid][7] = jmotif
								junction_lib[kid][8] = 1
							elif jmotif in set(jmotifs_rev):
								junction_lib[kid][6] = "-"
								junction_lib[kid][7] = jmotif
								junction_lib[kid][8] = 1
							else:
								junction_lib[kid][6] = "."
								junction_lib[kid][7] = jmotif
								junction_lib[kid][8] = 0
							
					if args.fix_bam:
						# use the motif orientation to add/replace the XS tag for this
						# alignment and write it out to the output BAM/SAM file
						pass
						
					
					# increment count for this junction
					junction_lib[kid][3] += 1
					
					# set the left side tail
					junction_lib[kid][4] = max(cu_match_length, junction_lib[kid][4])
					
					if junctions_in_aln > 1:
						# set right side tail for last junction
						junction_lib[kid_last][5] = max(cu_match_length, junction_lib[kid_last][5])
					
					# reset cumulative match length
					cu_match_length = 0 
					kid_last = kid
				
				# increment position for the next bit if this is a match, deletion or intron
				if aln.cigar[i][0] == 0 or aln.cigar[i][0] == 2 or aln.cigar[i][0] == 3: 
					ppos += aln.cigar[i][1]
			
			# saminished with the alignment's cigar so update the right tail of the last junction 
			if junctions_in_aln > 0:
				junction_lib[kid_last][5] = max(cu_match_length, junction_lib[kid_last][5])
			
			if aln_idx % 1000000 == 0:
				sys.stderr.write("[main] %d alignments parsed\n" % aln_idx)
	
#		if len(junction_lib.keys()) > 10:
#			break
		
	samin.close()
	fin.close()
	
	if False:
		# parse the FAI file into memory
		fin = open(args.fasta + ".fai", "r")
		for szl in fin:
			ll = szl.strip().split("\t")		
			dfai[ll[0]] = map(int, ll[1:])
		
		fin.close()
		
		# open the fasta
		fin = open(args.fasta, "r")
	
	# open output files
	fout_info = open(jinfo, "w")
	fout_bed = open(jbed, "w")
	
	# write header to the bed file (for the genome browser)
	fout_bed.write("track name=\"%s\" description=\"%s junctions\"\n" % (args.alignments, args.alignments))
	
	# -- look up the intron motifs from the discovered junctions
	sys.stderr.write("[main] printing junctions output to files\n")
	
	i = 0
	for kid in sorted(junction_lib.keys()):		
		jl = junction_lib[kid]
		
		jid = "JUNC{:08d}".format(i)
		i += 1

		fout_info.write("\t".join(map(str, jl)) + "\n")
		fout_bed.write(get_junc_bed_line(jid, jl) + "\n")
		
		if False:
			tid = jl[0]
			jstart = int(jl[1])-1
			jend = int(jl[2])
			jlen = jend-jstart
			
			if jlen < 4:
				# not long enough to have a motif...wtf, just print it back out
				jl += [".", "-", 0]
				fout_info.write("\t".join(map(str, jl)) + "\n")
				fout_bed.write(get_junc_bed_line(jid, jl) + "\n")
				# print "\t".join(map(str, jl))
			else:
			
				if False:
					# gather fai info for the reference
					roffset = dfai[tid][1]
					lnr = dfai[tid][2]
					lne = dfai[tid][3]
					
					# find how many rows, based on the lnr value, to the start of this 
					# junction sequence
					line_start = (jstart)/lnr
					line_start_offset = jstart - line_start*lnr
					line_start_rem = lnr-line_start_offset
					
					# offset to start of junction will be the reference offset + the junction position + 
					# the number of lines * the number of extra characters at the end of each line
					joffset = roffset + line_start*lne + line_start_offset
					
					# find number of positions remaining in the line the junction starts
					# line_rem_pos = jstart - num_lines_actual*lnr
					# subtract these bases from the length of the junction and find the number of lines we 
					# need to read in to capture the full sequence
					num_lines_junc = (jlen - line_start_rem)/lnr
					final_rem_pos = jlen - line_start_rem - num_lines_junc*lnr
					
					# find total read length
					total_read_length = (num_lines_junc)*lne + final_rem_pos
					if line_start_rem > 0:
						total_read_length += (line_start_rem + (lne-lnr))
					
					fin.seek(joffset)
					seq = fin.read(total_read_length)
					lseq = seq.split("\n")
					seq = "".join(lseq)
					# seq = seq[0:jlen]
					
					if len(seq) != jlen:
						sys.stderr.write("[main] Error: extracted sequence length is not equal to target length\n")
						print len(seq), jlen
			#			print line_rem_pos, num_lines_junc, final_rem_pos
			#			fin.close()
			#			return 1
					
					lseq = list(seq.upper())
				
				lseq = fetch_seq(fin, dfai, tid, jstart, jend)
				
				try:
					motif = [lseq[0], lseq[1], lseq[-2], lseq[-1]]
				except IndexError:
					print lseq
					print "\t".join(map(str, jl))
					return 1
				motif = "".join(motif)
				motif_rev = revcomp(motif)
				
				if motif in set(jmotifs):
					jl.append("+")
					jl.append(motif)
					jl.append(1)
				elif motif_rev in set(jmotifs):
					jl.append("-")
					jl.append(motif_rev)
					jl.append(1)
				else:
					jl.append(".")
					jl.append(motif)
					jl.append(0)
	
				fout_info.write("\t".join(map(str, jl)) + "\n")
				fout_bed.write(get_junc_bed_line(jid, jl) + "\n")
			
#			print "\t".join(map(str, jl))

#	fin.close()
	fout_info.close()
	fout_bed.close()
	
	return 0

def fetch_seq(fa, fai, rid, tstart, tend):
	# -- 
	# fetched sequencing from reference rid starting at tstart and ending at 
	# tend. fai should be a hashed table of the fasta index associated with 
	# the fasta file (pointer = fa)
	# --
	
	lseq = ""
	tlen = tend-tstart
	
	# gather fai info for the reference
	roffset = fai[rid][1]
	lnr = fai[rid][2]
	lne = fai[rid][3]
	
	# find how many rows, based on the lnr value, to the start of this 
	# junction sequence
	line_start = (tstart)/lnr
	line_start_offset = tstart - line_start*lnr
	line_start_rem = lnr-line_start_offset
	
	# offset to start of junction will be the reference offset + the junction position + 
	# the number of lines * the number of extra characters at the end of each line
	joffset = roffset + line_start*lne + line_start_offset
	
	# find number of positions remaining in the line the junction starts
	# line_rem_pos = tstart - num_lines_actual*lnr
	# subtract these bases from the length of the junction and find the number of lines we 
	# need to read in to capture the full sequence
	num_lines_junc = (tlen - line_start_rem)/lnr
	final_rem_pos = tlen - line_start_rem - num_lines_junc*lnr
	
	# find total read length
	total_read_length = (num_lines_junc)*lne + final_rem_pos
	if line_start_rem > 0:
		total_read_length += (line_start_rem + (lne-lnr))
	
	fa.seek(joffset)
	seq = fa.read(total_read_length)
	lseq = seq.split("\n")
	seq = "".join(lseq)
	# seq = seq[0:tlen]
	
	if len(seq) != tlen:
		sys.stderr.write("[main] Error: extracted sequence length is not equal to target length\n")
		print len(seq), tlen
	#			print line_rem_pos, num_lines_junc, final_rem_pos
	#			fin.close()
	#			return 1
	
	return seq
	

def get_junc_bed_line(jid, jl):
	
	jbed = []
	
	jbed.append(jl[0])
	jbed.append(jl[1]-jl[4]-1)
	jbed.append(jl[2]+jl[5])
	jbed.append(jid)
	jbed.append(jl[3])
	jbed.append(jl[6])
	jbed += [jbed[1], jbed[2], "255,0,0", 2]
	jbed.append("%d,%d" % (jl[4], jl[5]))
	jbed.append("0,%d" % (jl[2]-jl[1]+jl[4]+1))
	
	return "\t".join(map(str,jbed))

#
# reverse compliment a sequence
#
def revcomp(seq):
	
	lseq = list(seq.upper())
	i = 0
	j = len(lseq)-1
	
	tbl = { "A":"T", "C":"G", "G":"C", "T":"A", "N":"N" }
	
	while j >= i:
		tmp = tbl[lseq[i]]
		lseq[i] = tbl[lseq[j]]
		lseq[j] = tmp
		j -= 1
		i += 1
	
	return "".join(lseq) 


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


parser = argparse.ArgumentParser(description="Builds a junction library from BAM/SAM alignments.")
parser.add_argument('alignments', type=str, help="BAM alignments (or SAM with -S)")
parser.add_argument("fasta", type=str, help="FASTA sequence for the reference the alignments are relative to")
parser.add_argument("-S", dest="S", action="store_const", const=True, default=False, 
				help="Alignments are SAM (default expected: BAM)")
parser.add_argument("-o", dest="o", action="store", default="junctions", 
				help="Stub for output files (default: junctions)")
parser.add_argument("--fix-bam", dest="fix_bam", action="store_const", const=True, default=False, 
				help="Add the XS attribute to spliced alignments indicating the orientation of the splice.")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
