#!/usr/bin/env python
#
# tbam-to-gbam.py
#
# Shawn Driscoll
# 20121206
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script translates a BAM file of transcriptome alignments (currently only supporting
# ungapped alignments) into genomic alignments.
#

import sys
import subprocess as sp
import re
import argparse

def main(args):
	
	#
	# main script
	#

	# grab file names
	gtf = args.gtf
	bam = args.bam
	fasta = args.genome_fa
	
	bam_out = bam.split(".")[0] + ".genome"
	
	#
	# check files
	#
	try:
		fin = open(gtf,'r')
	except IOError,e:
		sys.stderr.write("Error: unable to open GTF file ({:s}). Bailing out.\n".format(gtf))
		return 1

	fin.close()

	try:
		fin = open(bam,'r')
	except IOError,e:
		sys.stderr.write("Error: unable to open BAM file ({:s}). Bailing out.\n".format(bam))
		return 1

	fin.close()
	
	try:
		fin = open(fasta,'r')
	except IOError,e:
		sys.stderr.write("Error: unable to open genome reference file ({:s}). Bailing out.\n".format(fasta))
		return 1

	fin.close()
	
	#
	# move on
	#
	
	#
	# build SAM header for genome
	#
	
	
	# first check for a fasta index in the same folder as the reference fasta file
	sys.stderr.write("> building SAM header for reference\n")
	try:
		fin = open(fasta + ".fai", 'r')
	except IOError,e:
		sys.stderr.write("> warning: there is no index for the FASTA specified. Samtools will be used to build one.\n")
		p1 = sp.Popen("samtools faidx {:s}".format(fasta).split())
		p1.wait()
		fin = open(fasta + ".fai", 'r')
	
	# parse the fasta index and build the header
	sam_header = []
	sam_header.append("@HD\tVN:1.0\tSO:unsorted")
	for szl in fin:
		szl = szl.strip()
		if len(szl) > 0:
			arl = szl.split("\t")
			sam_header.append("@SQ\tSN:" + arl[0] + "\tLN:" + arl[1])
	fin.close()
	
	# append program
	sam_header.append("@PG\tID:tbam-to-gbam.py\tVN:1.0\tCL:dunno")
	# append comment
	sam_header.append("@CO\tAlignments translanted from transcriptome alignments")
		
	#
	# parse the gtf
	#
	
	thash_genome = {}
	thash_lengths = {}
	thash_ref = {}
	
	# open the GTF file
	fin = open(gtf,"r")
	
	# parse GTF by transcript ID. also record total transcript lengths.
	sys.stderr.write("> parsing GTF\n")
	for szl in fin:
		szl = szl.strip()
		if len(szl) > 0:
			# line has content
			arl = szl.split("\t")
			
			# parse attributes field
			attrs = parse_gtf_attr(arl[8])
			tid = attrs['transcript_id']
	
			# add to hash
			if tid not in thash_genome:
				thash_genome[tid] = []
				thash_lengths[tid] = 0
			#<
	
			# insert exon
			thash_genome[tid].append([arl[0], int(arl[3]), int(arl[4]), arl[6]])
	
			# make sure the exons are sorted by position
			i = len(thash_genome[tid])
			while i > 1:
				i -= 1
				if thash_genome[tid][i][1] < thash_genome[tid][i-1][1]:
					# swap
					temp = list(thash_genome[tid][i-1])
					thash_genome[tid][i-1] = list(thash_genome[tid][i])
					thash_genome[tid][i] = list(temp)
				#<
			#<
	
			thash_lengths[attrs['transcript_id']] += int(arl[4])-int(arl[3])+1
		#<
	#<
	
	fin.close()
	
	#
	# create a second transcript hash with positions in terms of the reference
	#
	for tid in thash_genome.keys():
	
		# get set of exons for this transcript
		eset = thash_genome[tid]
		thash_ref[tid] = []
		ppos = 1
		for i in range(len(eset)):
			# create new exon coordinates
			ex_new = [eset[i][0], ppos, ppos + eset[i][2]-eset[i][1]]
			thash_ref[tid].append(ex_new)
			# increment position for the next exon
			ppos = ex_new[2]+1
		#<
	#<
	
	
	#sys.exit(1)
	
	# read through BAM file and write new alignments out a new BAM file
	
	sys.stderr.write("> writing BAM file...\n")
	
	p1 = sp.Popen("samtools view {:s}".format(bam).split(),stdout=sp.PIPE)
	# write to temporary bam file
	p2 = sp.Popen("samtools view -bS -o temp.bam -".split(),stdin=sp.PIPE)
	
	
	aln_rev = False
	
	# write header to new bam file
	for szl in sam_header:
		p2.stdin.write(szl + "\n")
	
	# loop through bam and translate alignments
	for szl in p1.stdout:
		szl = szl.strip()
		if len(szl) > 0:
			arl = szl.split("\t")
	
			tid = arl[2]
			ppos_ref = int(arl[3])
			ref_cigar = arl[5]
			types,lengths = parse_cigar(ref_cigar)
	
			# what is the total aligned length?
			aligned_length = 0
			for i in range(len(types)):
				if types[i] == "M":
					aligned_length += lengths[i]
	
			# what is the strand of the transcript id hit by this alignment?
			if thash_genome[tid][0][3] == '-':
				#sys.stderr.write("revising alignment position from {:d}".format(ppos_ref))
				ppos_ref = thash_lengths[tid] - (ppos_ref + aligned_length) + 2
				#sys.stderr.write(" to {:d}\n".format(ppos_ref))
				# reverse compliment the read
				arl[9] = rev_comp(arl[9])
	
			if len(types) == 1:
				if types[0] == "M":
					# one continuous alignment
					aln_start_ref = ppos_ref
					aln_end_ref = ppos_ref + lengths[0] - 1
	
					i = 0
					while i < len(thash_ref[tid]):
						exon = thash_ref[tid][i]
						exon_genome = thash_genome[tid][i]
	
						if aln_end_ref > exon[1] and aln_start_ref < exon[2]:
							# overlap
							if aln_end_ref > exon[2]:
								# junction hit, start cigar off by calling the amount
								# of match in this exon
								length_1 = exon[2] - aln_start_ref + 1
								cigar_new = str(length_1) + "M"
								# append the gap between exons
								gap = thash_genome[tid][i+1][1] - exon_genome[2] - 1
								if gap < 0:
									#print "Error: negative gap!"
									#print szl
									#print exon_genome,thash_genome[tid][i+1]
									sys.stderr.write("Error: negative gap!\n{:s}\n".format(szl))
									return 1
	
								cigar_new += str(thash_genome[tid][i+1][1] - exon_genome[2] - 1) + "N"
	
								# is the rest of the alignment contained in this exon?
								if aln_end_ref <= thash_ref[tid][i+1][2]:
									# yes so we can just report the remaining length here
									length_2 = lengths[0] - length_1
									cigar_new += str(length_2) + "M"
	
									arl[5] = cigar_new
									arl[3] = ppos_ref - exon[1] + exon_genome[1]
									arl[2] = exon_genome[0]
									#print "\t".join(map(str,arl))
									p2.stdin.write("\t".join(map(str,arl)) + "\n")
	
								else:
									# no, this alignment segment will be the length of this exon
									length_2 = thash_ref[tid][i+1][2]-thash_ref[tid][i+1][1]+1
									cigar_new += str(length_2) + "M"
									sys.stderr.write("> Warning: this read hits more than 2 exons. {:d} bases remaining\n".format(lengths[0]-length_1-length_2))
	
	
							else:
								# alignment is contained within this exon
								arl[5] = str(lengths[0]) + "M"
								arl[3] = ppos_ref - exon[1] + exon_genome[1]
								arl[2] = exon_genome[0]
								#print "\t".join(map(str,arl))
								p2.stdin.write("\t".join(map(str, arl)) + "\n")
	
							i = len(thash_ref[tid])
	
						i += 1
	
			else:
				sys.stderr.write("> Warning: dropped a read with cigar {:s}\n".format(ref_cigar))
	
	p2.stdin.close()

	#
	# TODO: this script should include some way to parse the genome alignments to remove duplicates. in many cases the
	# same read with the same alignment will be repeated in this output thanks to genomes with alternative splicing.
	# I tried passing the SAM data through 'uniq' via a pipe but that didn't fix all of the duplicates. 
	#	

	# 
	# it looks like this could be accomplished by name-sorting the alignments and passing them through 'uniq'.
	# this will take hell time I'm sure
	#
	sys.stderr.write("> sorting BAM by read name\n")
	p1 = sp.Popen("samtools sort -n temp.bam temp-n".split())
	p1.wait()
	
	sys.stderr.write("> removing redundant alignments\n")
	
	# open pipe to convert output to BAM format
	p2 = sp.Popen("samtools view -bS -o temp.bam -".split(),stdin=sp.PIPE)
	
	# pass header to new file
	p1 = sp.Popen("samtools view -H temp-n.bam".split(),stdout=sp.PIPE)
	
	for szl in p1.stdout:
		szl = szl.strip()
		if len(szl) > 0:
			p2.stdin.write(szl + "\n")
	
	# pass the name sorted data through uniq
	p1 = sp.Popen("samtools view temp-n.bam".split(), stdout=sp.PIPE)
	p3 = sp.Popen("uniq".split(), stdin=p1.stdout, stdout=sp.PIPE)
	
	for szl in p3.stdout:
		szl = szl.strip()
		if len(szl) > 0:
			p2.stdin.write(szl + "\n")
	
	p2.stdin.close()
	
	
	# bam file should be complete now we can sort it if requested
	if args.no_sort:
		sp.call("mv temp.bam {:s}.bam".format(bam_out),shell=True)
	else:
		sys.stderr.write("> sorting new bam file\n")
		sp.call("samtools sort temp.bam {:s}".format(bam_out),shell=True)
		sp.call("rm temp.bam", shell=True)
	
	sp.call("rm temp-n.bam", shell=True)
	
	sys.stderr.write("> done!\n")
	
#}

def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	#
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

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

#<


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

#<

#
# MAIN ENTRY POINT
#

parser = argparse.ArgumentParser(description="Translate transciptome alignments in BAM format to genomic coordinates in BAM format.")

parser.add_argument('genome_fa', type=str, action='store', help='Genome reference in FASTA.')
parser.add_argument('gtf', type=str, action='store', help='GTF file describing the transcriptome to which the data was aligned.')
parser.add_argument('bam', type=str, action='store', help='BAM file containing alignments to translate.')
parser.add_argument('--no-sort', dest='no_sort', action='store_const', const=True, default=False, help="Do not sort the final BAM file (default: sort)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
#}
			




