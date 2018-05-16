#!/usr/bin/env python
#==============================================================================
# make-junc-db.py
#
# Shawn Driscoll
# 20150102
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Make junction database (include exon/intron boundaries) for mapping to 
# splices and retained introns. 
#==============================================================================

import sys, argparse, re
import subprocess as sp
from os import unlink
from os.path import isfile, expanduser

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r


HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	juncdb = {}
	retdb = {}
	tid = ""
	gid = ""
	gname = ""
	ltid = ""
	laln = []
	jid = ""
	stub = ""
	wwin = args.l - args.o

	if wwin <= 0:
		sys.stderr.write("[main] error: read length minus overlap results in 0 or negative length\n")
		return 1

	# check input file
	if not file_exists(args.gtf):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.gtf)
		return 1

	# setup output file name stub

	# strip off the folder path
	tmp0 = (args.gtf).split("/")
	# split the file name (last element in the list)
	tmp = tmp0[-1].split(".")
	# stick the file name back together with all but the last dot
	stub = ".".join(tmp[0:len(tmp)-1])
	tmp = stub + "_psiTheta_{}".format(args.l)
	stub = tmp

	juncdb_gtf = stub + ".gtf"
	juncdb_fa = stub + ".fa"
	juncdb_index = stub + ".index.fa"


	#
	# Start parsing the GTF
	#

	sys.stderr.write("Parsing {}\n".format(args.gtf))
	fin = open(args.gtf, "r")
	for szl in fin:
		aln = szl.strip().split("\t")

		# skip any non exon lines
		if aln[2] != "exon":
			continue

		# find the transcript id
		rres = re.search("transcript_id \"([^\"]+)\"", aln[8])
		if rres:
			tid = rres.group(1)
		else:
			sys.stderr.write("[main] failed to parse transcript id")
			fin.close()
			return(1)

		# check if current transcript id matches the last one. if so then we have a junction
		# because we have gone from one exon to the next
		if tid == ltid:
			# junction is from laln[4] to aln[3]
			jid = "{}:{}-{}".format(aln[0], str(int(laln[4])+1), str(int(aln[3])-1))
			# is this junction in the database already?
			if jid not in juncdb:
				juncdb[jid] = dict(ref=aln[0], strand=aln[6], j1=[max(int(laln[3]), int(laln[4])-wwin+1), int(laln[4])], 
					j2=[int(aln[3]), min(int(aln[4]), int(aln[3])+wwin-1)], gnames={}, gids={}, tids={})

			# add transcript id to the set
			juncdb[jid]["tids"][tid] = 0

			# add gene id to the set
			rres = re.search("gene_id \"([^\"]+)\"", aln[8])
			gid = ""
			if rres:
				gid = rres.group(1)
				juncdb[jid]['gids'][rres.group(1)] = 0

			# add gene name to the set
			rres = re.search("gene_name \"([^\"]+)\"", aln[8])
			gname = ""
			if rres:
				gname = rres.group(1)
				juncdb[jid]['gnames'][rres.group(1)] = 0

			#
			# exon edges from current and last will be "retention" ids
			#
			
			# do previous (5p or donor)

			rid = "{}:{}:5p".format(aln[0], str(int(laln[4])+1))
			if rid not in retdb:
				retdb[rid] = dict(ref=aln[0], strand=aln[6], 
					rng=[max(int(laln[3]), int(laln[4])-wwin+1), min(int(laln[4])+wwin, int(aln[4]))], 
					gnames={}, gids={}, tids={})

			retdb[rid]['tids'][tid] = 0
			if len(gid) > 0:
				retdb[rid]['gids'][gid] = 0
			if len(gname) > 0:
				retdb[rid]['gnames'][gname] = 0

			# do current (3p or acceptor)

			rid = "{}:{}:3p".format(aln[0], str(int(aln[3])-1))
			if rid not in retdb:
				retdb[rid] = dict(ref=aln[0], strand=aln[6], 
					rng=[max(int(laln[3]), int(aln[3])-wwin), min(int(aln[3])+wwin-1, int(aln[4]))], 
					gnames={}, gids={}, tids={})

			retdb[rid]['tids'][tid] = 0
			if len(gid) > 0:
				retdb[rid]['gids'][gid] = 0
			if len(gname) > 0:
				retdb[rid]['gnames'][gname] = 0


		laln = list(aln)
		ltid = tid

	fin.close()

	# last line
	if tid == ltid:
		# junction is from laln[4] to aln[3]
		jid = "{}:{}-{}".format(aln[0], laln[4], aln[3])
		if jid not in juncdb:
			juncdb[jid] = dict(ref=aln[0], strand=aln[6], j1=[max(int(laln[3]), int(laln[4])-wwin+1), int(laln[4])], 
				j2=[int(aln[3]), min(int(aln[4]), int(aln[3])+wwin-1)], gnames={}, gids={}, tids={})

		juncdb[jid]["tids"][tid] = 0

		rres = re.search("gene_id \"([^\"]+)\"", aln[8])
		gid = ""
		if rres:
			gid = rres.group(1)
			juncdb[jid]['gids'][rres.group(1)] = 0

		rres = re.search("gene_name \"([^\"]+)\"", aln[8])
		gname = ""
		if rres:
			gname = rres.group(1)
			juncdb[jid]['gnames'][rres.group(1)] = 0

		# exon edges from current and last will be "retention" ids
		
		# do previous (5p or donor)

		rid = "{}:{}:5p".format(aln[0], str(int(laln[4])+1))
		if rid not in retdb:
			retdb[rid] = dict(ref=aln[0], strand=aln[6], 
				rng=[max(int(laln[3]), int(laln[4])-wwin+1), min(int(laln[4])+wwin, int(aln[4]))], 
				gnames={}, gids={}, tids={})

		retdb[rid]['tids'][tid] = 0
		if len(gid) > 0:
			retdb[rid]['gids'][gid] = 0
		if len(gname) > 0:
			retdb[rid]['gnames'][gname] = 0

		# do current (3p or acceptor)

		rid = "{}:{}:3p".format(aln[0], str(int(aln[3])-1))
		if rid not in retdb:
			retdb[rid] = dict(ref=aln[0], strand=aln[6], 
				rng=[max(int(laln[3]), int(aln[3])-wwin), min(int(aln[3])+wwin-1, int(aln[4]))], 
				gnames={}, gids={}, tids={})

		retdb[rid]['tids'][tid] = 0
		if len(gid) > 0:
			retdb[rid]['gids'][gid] = 0
		if len(gname) > 0:
			retdb[rid]['gnames'][gname] = 0

	sys.stderr.write("writing {}...\n".format(juncdb_gtf))
	fout = open(juncdb_gtf, "w")

	# write out as a GTF where each "transcript" is a single junction event or pair of features
	for jid in sorted(juncdb.keys()):
		oId = ",".join(sorted(juncdb[jid]['tids'].keys()))
		tid = jid
		gname = ",".join(sorted(juncdb[jid]['gnames'].keys()))
		gid = ",".join(sorted(juncdb[jid]['gids'].keys()))

		# first line
		lout = [juncdb[jid]['ref'], "psi", "exon", 
			str(juncdb[jid]['j1'][0]), str(juncdb[jid]['j1'][1]), ".", juncdb[jid]['strand'], ".", 
			"gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; oId \"{}\";".format(gid, tid, gname, oId)]

		fout.write("\t".join(lout) + "\n")

		# second line
		lout = [juncdb[jid]['ref'], "psi", "exon", 
			str(juncdb[jid]['j2'][0]), str(juncdb[jid]['j2'][1]), ".", juncdb[jid]['strand'], ".", 
			"gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; oId \"{}\";".format(gid, tid, gname, oId)]

		fout.write("\t".join(lout) + "\n")

	for rid in sorted(retdb.keys()):
		oId = ",".join(sorted(retdb[rid]['tids'].keys()))
		tid = rid
		gname = ",".join(sorted(retdb[rid]['gnames'].keys()))
		gid = ",".join(sorted(retdb[rid]['gids'].keys()))

		# only one line
		lout = [retdb[rid]['ref'], "theta", "exon", 
			str(retdb[rid]['rng'][0]), str(retdb[rid]['rng'][1]), ".", retdb[rid]['strand'], ".", 
			"gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; oId \"{}\";".format(gid, tid, gname, oId)]

		fout.write("\t".join(lout) + "\n")

	sys.stderr.write("done\n")

	fout.close()


	# use gffread to extract the sequences from the supplied genome file
	cmd = "gffread -w {} -g {} {}".format(stub + ".fa", args.fa, stub + ".gtf")
	runcmd(cmd)

	if not args.basic_index:

		# kmer the reference
		cmd = "kmercountexact.sh in={} k={} tuc=t out=temp_mers.fa overwrite=t".format(juncdb_fa, args.l)
		runcmd(cmd)
		
		# match back to the reference with seal
		cmd = "seal.sh in=temp_mers.fa ref={} k={} rcomp=t rename=t tuc=t hdist=0 mm=f trd=t out=temp_mers_matched.fa overwrite=t".format(juncdb_fa, args.l)
		runcmd(cmd)
		
		# now run that "other" python script to make the final index
		cmd = "{}/coding/python/juncdb-filter-matched-ref-mers.py {} {}".format(HOME, juncdb_gtf, "temp_mers_matched.fa")
		# open output file
		fout = open(juncdb_index, "w")
		# run the command 
		sys.stderr.write("CMD: {}\n".format(cmd))
		p1 = sp.Popen(cmd.split(), stdout=fout)
		p1.wait()
		fout.close()
		
		# clean up
		unlink("temp_mers.fa")
		unlink("temp_mers_matched.fa")
		
		# compress 
		runcmd("gzip {}".format(juncdb_fa))
		runcmd("gzip {}".format(juncdb_index))

#	# make fai
#	cmd = "samtools faidx {}".format(stub + ".fa")
#	sys.stderr.write("CMD: {}\n".format(cmd))
#	pid = sp.Popen(cmd.split())
#	pid.wait()

	if args.bowtie:
		# make bowtie index
		cmd = "bowtie-build -o 1 {} {}".format(stub + ".fa", stub)
		ferr = open("/dev/null", "w")
		sys.stderr.write("CMD: {}\n".format(cmd))
		pid = sp.Popen(cmd.split(), stderr=ferr, stdout=ferr)
		pid.wait()
		ferr.close()

	return 0


# --
# runcmd
# run a system level command in subprocess. optionally you can return the process.
# if the process isn't returned then the function waits for the process to finish
def runcmd(cmd, returnProcess=False):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)

def file_exists(fname):
	try:
		fin = open(fname)
		fin.close()
	except IOError as e:
		return False

	return True

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="...")
parser.add_argument('gtf', type=str, help="GTF to extract junctions from")
parser.add_argument('fa', type=str, help="Genome FASTA reference")
parser.add_argument('-o', type=int, action="store", default=4, 
	help="Base pairs of overlap you'll want to require for valid hits [4]")
parser.add_argument('-l', type=int, action="store", default=50, 
	help="Read length [50]")
parser.add_argument('--bowtie', action="store_const", const=True, default=False, 
	help="Build bowtie index [off]")
parser.add_argument('--basic-index', action="store_const", const=True, default=False, 
	help="Basic index is only the first FASTA sequence file with the basic junction and exon-intron boundary sequences.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

