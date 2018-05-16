#!/usr/bin/python

# run seal to intersect a reference with itself

import os, sys, argparse, re
import subprocess as sp
from os.path import isfile
from math import floor

def main(args):


	g = {} # hash that will contain paired sequences and their relative kmer-similarity
	g0 = {}
	kcontrol = {} # hash containing the number of kmers per sequence
	# output file name
	outname = "{}.mapp".format(args.fa)

	dresult = {}
	dwritten = {}
	np_kept = 0
	np_found = 0

	# fetch the length of the sequences from the fasta's FAI index
	dtlens = {}
	if not isfile("{}.fai".format(args.fa)):
		sys.stderr.write("We need to create an fai index for the input FASTA\n")
		cmd = "samtools faidx {}".format(args.fa)
		p1 = sp.Popen(cmd.split())
		p1.wait()

	sys.stderr.write("Reading {}.fai\n".format(args.fa))
	fin = open("{}.fai".format(args.fa), "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		dtlens[aln[0]] = float(aln[1])
		dresult[aln[0]] = [aln[0], "none", 0, 0, float(aln[1]), args.k, floor((float(aln[1])-args.k+1)/args.s)]

	fin.close()

	gpairs = []

	# build seal command
	seal_cmd = "seal.sh in={} ref={} k={} trd=t mm=f rename=t rcomp=f ambig=all out=stdout.fa qskip={}".format(args.fa, args.fa, args.k, args.s)
	if args.d > 0:
		seal_cmd = seal_cmd + " qhdist={}".format(args.d)
	if args.p > 0:
		seal_cmd = seal_cmd + " threads={}".format(args.p)

	# start subprocess with seal
	ferr = open("/dev/null", "w")
	p1 = sp.Popen(seal_cmd.split(), stdout=sp.PIPE, stderr=ferr)

	sys.stderr.write("CMD: " + seal_cmd + "\n")

	# read in results
	for szl in p1.stdout:
		if szl[0] == ">":
			# only parse reference name lines
			aln = szl.strip().split("\t")
			# base sequence is the first one
			base = re.sub("^>", "", aln[0])
			# start a list in the hash for this reference
			g[base] = []
			# note each sequence this one has kmer matches to and how many
			for i in range(1, len(aln)):
				tmp = aln[i].split("=")

				if tmp[0] == base:
					# this indicates number of hits to itself, keep this for later when we calculate ratios
					# and leave it out of the targets list for this reference
					kcontrol[base] = float(tmp[1])
				else:
					# note the sequence and number of kmer matches for the current base to the current
					# target
					g[base].append([tmp[0], float(tmp[1]), 0])


	p1.stdout.close()
	ferr.close()

	sys.stderr.write("Matched {} sequences (if some are missing it may be due to length relative to k)\n".format(len(g.keys())))

	sys.stderr.write("Generating pairwise matches.\n")

	# loop back through all of the references. if they have targets then calculate ratios

	for tid in g.keys():
		nhits = len(g[tid])
		if nhits > 0:

			# loop through and generate the pairs
			for lg in g[tid]:
				np_found += 1

				krat = lg[1]*1.0/kcontrol[tid]
				if krat > args.t:
					np_kept += 1
					lout = list(dresult[tid])
					lout[1] = lg[0]
					lout[2] = krat
					lout[3] = lg[1]
					gpairs.append(lout)

				if np_found % 100000 == 0:
					sys.stderr.write("found {} pairings. {} over threshold\n".format(np_found, np_kept))

#		else:
#			# enter a row for the transcript so we see it in the output
#			gpairs.append(list(dresult[tid]))

	sys.stderr.write("found {} pairings. {} over threshold\n".format(np_found, np_kept))
	sys.stderr.write("Writing report to {}\n".format(outname))

	# print header
	fout = open(outname, "w")
	fout.write("reference\ttarget\tshared_ratio\tshared_k_num\treference_length\tk\treference_num_k\n")

	# write result
	for i in range(len(gpairs)):
		dwritten[gpairs[i][0]] = 0
		fout.write("\t".join(map(str, gpairs[i])) + "\n")

	# write out any that didn't get handled.
	for tid in dresult.keys():
		if tid not in dwritten:
			fout.write("\t".join(map(str, dresult[tid])) + "\n")

	fout.close()

	return 0





#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(
	description="""...""")

parser.add_argument("fa", type=str, 
	help="FASTA reference")

# quantification options
parser.add_argument("-k", type=int, dest="k", default=31, action="store", 
	help="mer length [31]")
parser.add_argument("-s", type=int, dest="s", default=2, action="store", 
	help="skip this many query kmers [2]")
parser.add_argument("-p", type=int, dest="p", default=0, action="store",
	help="number of threads for seal [all]")
parser.add_argument("-d", type=int, dest="d", default=0, action="store",
	help="allowed mismatches when matching kmers [0]")
parser.add_argument('-t', type=float, default=0, 
	help="Threshold for similarity. Only connections above this threshold are retained. Enter 0 to keep all. [0]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

