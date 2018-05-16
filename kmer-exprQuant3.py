#!/usr/bin/env python
#==============================================================================
# kmer-exprQuant2.py
#
# Shawn Driscoll
# 20160128
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This version counts reads/fragments instead of only mers and considers
# best hits, discarding those with fewer kmer matches than best. ties are 
# split according to the typical rules.
#==============================================================================

import sys, argparse, math, re, os
import subprocess as sp
from os.path import isfile
import hashlib

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
import numpy as np
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
	nassigned = 0
	current_read = ""
	last_read = ""
	logfile = "keq.log"

	# -------------------------------------------------------------------------
	# load GTF annotation
	# -------------------------------------------------------------------------

	sys.stderr.write("loading GTF annotation...\n")
	annot = parse_gtf(args.gtf)
	sys.stderr.write("found {} features\n".format(len(annot.keys())))

	# -------------------------------------------------------------------------
	# prepare reads for seal
	# -------------------------------------------------------------------------

	# if PE then merge the files into a single interleaved file. if stranded then 
	# flip the mate

	pe = False
	if len(args.reads) > 1:
		pe = True

	istranded = 0
	if args.r_stranded or args.f_stranded:
		istranded = 1

	if pe:
		cmd = "reformat.sh in={} in2={} out={}_reads.fa.gz pigz=t".format(args.reads[0], args.reads[1], stub)

		if istranded != 0:
			cmd = cmd + " rcompmate=t"

		# run it
		rres = runcmd(cmd, False)
		fseal_in = "{}_reads.fa.gz".format(stub)
	else:
		fseal_in = args.reads[0]

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

	# write seal output to log
	logout = open(logfile, "w")

	cmd = "seal.sh in={} ref={} k={} rename=t trd=t ambig=all hdist=0 mm=t overwrite outm={}".format(fseal_in, refin, args.k, fifo_seal)
	if pe:
		cmd += " kpt=t interleaved=t"
	if istranded != 0:
		cmd += " rcomp=f"

	cmd += " qskip={}".format(args.kmer_interval)

	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split(), stderr=logout)

	fin = open(fifo_seal, "r")

	for szl in fin:

		# PARSE SEAL OUTPUT

		rres = re.search("^>", szl)
		if rres:
			nlines += 1

			# name row, parse this for hits
			szl = re.sub('^>', '', szl)
			aln = szl.strip().split("\t")

			current_read = aln[0]
			if pe and current_read == last_read:
				# skip second mates because info is redundant
				continue

			# track gene ids hit
			tgid = {}

			# expand the targets
			targets = []
			targets0 = []
			khits = []

			# apply some more thinking - if a read hits multiple features but has a single "best"
			# then that's what we'll call the actual target. if the read has multiple "best" then
			# the alignment will be shared.
			for i in range(1, len(aln)):
				# the name consists of the 'tid' key for the annot dict as well as a hit count
				# that follows an '=' sign. so we can split on that

				# parse target name and hit count
				atmp = aln[i].split("=")
				atmp[1] = float(atmp[1])

				if atmp[1] >= args.min_kmers:
					targets0.append(list(atmp))
					khits.append(atmp[1])

			if len(targets0) > 0:

				if len(targets0)==1:
					# nothing to consider...
					tgid[annot[targets0[0][0]]['gid']] = 0
					targets.append(list(targets0[0]))

				else:

					# sort the hits
					khits_order = np.argsort(khits)

					i = len(khits_order)-1

					# keep the best one for sure...
					try:
						tgid[annot[targets0[khits_order[i]][0]]['gid']] = 0
						targets.append(list(targets0[khits_order[i]]))
					except IndexError as e:
						sys.stderr.write("{}\n".format(", ".join(map(str, khits_order))))
						sys.stderr.write("{}, {}, {}\n".format(i, len(khits_order), len(targets0)))
						#p1.terminate()
						return 1

					while i > 0:
						i -= 1
						if targets0[khits_order[i]][1]==targets0[khits_order[i+1]][1]:
							# same hit count, keep it as well...
							tgid[annot[targets0[khits_order[i]][0]]['gid']] = 0
							targets.append(list(targets0[khits_order[i]]))
						else:
							# lower hit count, time to drop out
							i = -1

				if len(targets) > 0:

					nassigned += 1
					
					# down-weighting for multiple targets and multiple gene ids
					tgid_weight = float(1)/len(targets) * float(1)/(len(tgid.keys())**2)

					# loop through targets and assign hits
					for i in range(len(targets)):
						hits = (1*tgid_weight)
						# first assign the hit to the annot table
						annot[targets[i][0]]['hits'] += hits
						# track total depth
						total_depth += hits

			if (nlines % 1000000) == 0:
				sys.stderr.write("matched / assigned / %: {} / {} / {:0.2f}\n".format(nlines, nassigned, float(nassigned)/nlines*100.0))

			# track read names
			last_read = current_read


	sys.stderr.write("matched / assigned / %: {} / {} / {:0.2f}\n".format(nlines, nassigned, float(nassigned)/nlines*100.0))

	sys.stderr.write("done. writing output to stdout...\n")

	fin.close()
	logout.close()

	#
	# remove temporary files
	#

	if isfile("{}_ref.fa".format(stub)):
		os.unlink("{}_ref.fa".format(stub))

	try:
		os.unlink(fifo_seal)
		os.unlink("{}_reads.fa.gz".format(stub))
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

parser = argparse.ArgumentParser(description="Quantfies reads against a reference using kmers.")
parser.add_argument('ref', type=str, help="FASTA reference")
parser.add_argument('gtf', type=str, help="GTF annotation corresponding to the reference")
parser.add_argument('reads', type=str, nargs='+', metavar='reads', 
	help="First or First+Second mate FASTQ files")

parser.add_argument('-k', action="store", type=int, default=31, 
	help="kmer length [31]")
parser.add_argument('-i', '--kmer-interval', action="store", default=2, type=int, 
	help="Query kmer interval (1 means all kmers are used) [2]")
parser.add_argument('-m', '--min-kmers', action="store", type=int, default=2, 
	help="Minimum number of kmer hits to a reference to count as a hit for that read [2]")
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

