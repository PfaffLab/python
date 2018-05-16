#!/usr/bin/env python
#==============================================================================
# template.py
#
# Shawn Driscoll
# 20151202
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Pipeline to count mer hits to psi and theta splice site sequences. "alignments"
# come from Seal and input will be a fasta file containing the targets hit. 
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib

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
	total_depth = 0
	stub = hashlib.md5(args.reads[0]).hexdigest()
	fifo_seal = "{}_sealout.fa".format(stub)
	nlines = 0

	# -------------------------------------------------------------------------
	# load GTF annotation
	# -------------------------------------------------------------------------

	sys.stderr.write("loading GTF annotation...\n")
	annot = parse_gtf(args.gtf)
	sys.stderr.write("found {} features\n".format(len(annot.keys())))

	# -------------------------------------------------------------------------
	# count read kmers
	# -------------------------------------------------------------------------

	# paired end
	istranded = 0
	if args.r_stranded or args.f_stranded:
		istranded = 1

	if len(args.reads) > 1:
		fseal_in = kmer_reads(args.reads[0], args.reads[1], args.k, stub, istranded)
	else:
		fseal_in = kmer_reads(args.reads[0], "", args.k, stub, istranded)

	# -------------------------------------------------------------------------
	# reverse complement the reference if reverse stranded
	# -------------------------------------------------------------------------

	refin = args.ref
	if args.r_stranded:
		cmd = "reformat.sh in={} rcomp out={}_ref.fa".format(args.ref, stub)
		rres = runcmd(cmd, False)
		refin = "{}_ref.fa".format(stub)

	# -------------------------------------------------------------------------
	# run Seal and parse output
	# -------------------------------------------------------------------------

	try:
		os.unlink(fifo_seal)
	except:
		# don't do anything
		pass

	sys.stderr.write("+mkfifo {}\n".format(fifo_seal))
	os.mkfifo(fifo_seal)

	cmd = "seal.sh in={} ref={} k={} rename=t trd=t ambig=all hdist=0 mm=t overwrite outm={}".format(fseal_in, refin, args.k, fifo_seal)
	if istranded != 0:
		cmd += " rcomp=f"
	sys.stderr.write("CMD: {}".format(cmd))
	p1 = sp.Popen(cmd.split())

	fin = open(fifo_seal, "r")

	for szl in fin:

		# PARSE SEAL OUTPUT

		rres = re.search("^>", szl)
		if rres:
			nlines += 1

			# name row, parse this for hits
			szl = re.sub('^>', '', szl)
			aln = szl.strip().split("\t")

			# track gene ids hit
			tgid = {}

			# expand the targets
			targets = []
			for i in range(1, len(aln)):
				# the name consists of the 'tid' key for the annot dict as well as a hit count
				# that follows an '=' sign. so we can split on that

				# parse target name and hit count
				atmp = aln[i].split("=")
				# name of the fasta row is the count because it was a kmer database
				atmp[1] = float(aln[0])

				# note gene id so we can figure out downweighting later
				tgid[annot[atmp[0]]['gid']] = 0

				# append to list of targets
				targets.append(list(atmp))

			# down-weighting for multiple targets and multiple gene ids
			tgid_weight = float(1)/len(targets) * float(1)/(len(tgid.keys())**2)

			# loop through targets and assign hits
			for i in range(len(targets)):
				hits = (targets[i][1]*tgid_weight)
				# first assign the hit to the annot table
				annot[targets[i][0]]['hits'] += hits
				# track total depth
				total_depth += hits

			if (nlines % 1000000) == 0:
				sys.stderr.write("parsed {} mer hits\n".format(nlines))


	sys.stderr.write("parsed {} mer hits\n".format(nlines))

	sys.stderr.write("done. writing output to stdout...\n")

	fin.close()

	#
	# remove temporary files
	#

	try:
		os.unlink(fifo_seal)
		os.unlink(fseal_in)
		os.unlink("{}_ref.fa".format(stub))
	except:
		pass

	#
	# hits are all assigned. now we can produce the output
	#

	# print header
	print "\t".join(["transcript_id", "gene_id", "gene_name", "chrom", "strand", "length", "khits"])

	for tid in sorted(annot.keys()):

		lout = [
			annot[tid]['tid'], 
			annot[tid]['gid'],
			annot[tid]['gname'],
			annot[tid]['chrom'],
			annot[tid]['strand'],
			str(annot[tid]['length']),
			"{:f}".format(annot[tid]['hits'])
		]

		print "\t".join(map(str, lout))

	return 0

# 
# this function deals with calling kmercountexact.sh on the input reads.
# in the reads are stranded or paired-end then different options are
# used.
def kmer_reads(r1, r2, k, stub, strand):
	
	# strand = 0 (unstranded), 1 (stranded)

	try:
		os.unlink("{}_mers.fa.gz".format(stub))
	except:
		pass

	pe = False
	cmd = ""
	kmercmd = "kmercountexact.sh k={} overwrite=t out={}_mers.fa.gz pigz=t".format(k, stub)
	infile = r1

	if not isfile(r1):
		sys.stderr.write("Cannot find first-mate reads file")
		return 1

	if len(r2) > 0:
		if not isfile(r2):
			sys.stderr.write("Cannot find second-mate reads file")
			return 1

		pe = True

	cmd = ""

	if strand != 0:
		kmercmd += " rcomp=f"

		if pe==True:
			# if PE is true and the data is stranded then the second mate reads 
			# have to be reverse complemented prior to kmering
			sys.stderr.write("[kmer_reads] reverse complementing second mates\n")
			cmd = "reformat.sh in={} in2={} out={}.fa.gz rcompmate overwrite pigz=t".format(r1, r2, stub)
			rres = runcmd(cmd, False)

			kmercmd += " in={}.fa.gz".format(stub)
		else:
			kmercmd += " in={}".format(r1)

	else:
		# count all kmers reverse and forward
		kmercmd += " in={}".format(r1)
		if pe==True:
			kmercmd += " in2={}".format(r2)

	sys.stderr.write("[kmer_reads] counting kmers...\n")
	rres = runcmd(kmercmd, False)

	# return the name of the file that will be input to seal
	return("{}_mers.fa.gz".format(stub))


#
# parser the GTF
def parse_gtf(fname):
	# variables
	gtfdb = {}
	szl = ""
	aln = []

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		if aln[2] != "exon":
			continue

		# parse attributes
		attr = parse_gtf_attr(aln[8])

		if attr['transcript_id'] not in gtfdb:
			gtfdb[attr['transcript_id']] = {'gid': attr['gene_id'], 'gname': attr['gene_name'], 
			'tid': attr['transcript_id'], 'hits': 0, 'strand': aln[6], 'chrom': aln[0], 'length': 0}
		
		# increment length
		gtfdb[attr['transcript_id']]['length'] += (int(aln[4])-int(aln[3]))+1

	fin.close()

	return gtfdb


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

def runcmd(cmd, returnProcess):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Quantfies reads against a reference using kmers")
parser.add_argument('ref', type=str, help="FASTA reference")
parser.add_argument('gtf', type=str, help="GTF annotation corresponding to the reference")
parser.add_argument('reads', type=str, nargs='+', metavar='reads', 
	help="First or First+Second mate FASTQ files")

parser.add_argument('-k', type=int, default=26, 
	help="kmer length [26]")
parser.add_argument('-r', '--r-stranded', action="store_const", const=True, default=False, 
	help="Reverse stranded library [off]")
parser.add_argument('-f', '--f-stranded', action="store_const", const=True, default=False, 
	help="Forward stranded library [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

