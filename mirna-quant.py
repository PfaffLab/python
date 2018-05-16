#!/usr/bin/python
#
# bwa-quant.py
#
# porting the perl version over to python...just because
#
#

import os,sys,re
import argparse
from hashlib import md5
import subprocess as sp

def main(args):

	pe = False
	hits = {}
	hfactor = 0
	stats = dict(aligned=0, assigned=0, unique=0, ambig=0, gene_ambig=0)
	ghits = {} # used for gene/locus hit ambiguity counting
	ttog = {} # for translating transcript id to gene/locus id for ambiguity
	use_ambig = False

	# seal command
	seal_base = "/home/colfax/opt/bbmap/seal.sh"

	# --
	# check options
	# --
	

	if len(args.B) > 0:
		try:
			fin = open(args.B, "r")
			fin.close()
		except:
			sys.stderr.write("Error: failed to open {} (-B)\n".format(args.B))
			return(1)

		use_ambig = True
		# open and parse the file
		fin = open(args.B, "r")
		for szl in fin:
			aln = szl.strip().split("\t")
			ttog[aln[0]] = aln[1]
#			hits[aln[0]] = [0, 0]
		fin.close()

	# -- 
	# load fai if it exists, otherwise generate it and load it to build the
	# initial hits hash so that all targets are known as well as their lengths
	# --

	sys.stderr.write("loading {}.fai to parse in all target names\n".format(args.ref))

	try:
		fin = open(args.ref + ".fai", "r")
	except IOError:
		# can't open it so make it
		p1 = sp.Popen("samtools faidx {}".format(args.ref).split())
		p1.wait()
		fin = open(args.ref + ".fai", "r")

	for szl in fin:
		aln = szl.strip().split("\t")
		hits[aln[0]] = [int(aln[1]), 0, 0, 0]


	# -- 
	# run seal to start kmer matching process to the reference
	# --
	
	# run alignment(s)
	if args.rand_sample > 1:
		cmd = "reformat.sh in={} srt={} out=rand_sample.fq t={} overwrite=t".format(args.reads, args.rand_sample, args.t)
		p1 = sp.Popen(cmd.split())
		p1.wait()
		cmd = seal_base + " in=rand_sample.fq overwrite=t out=stdout.fq forbidn=f ambiguous=all stats=seal.stats rename=t trd=t ref={} hdist={} k={} t={}".format(args.ref, args.d, args.k, args.t)
	elif args.rand_sample > 0 and args.rand_sample < 1:
		cmd = "reformat.sh in={} samplerate={} out=rand_sample.fq t={} overwrite=t".format(args.reads, args.rand_sample, args.t)
		p1 = sp.Popen(cmd.split())
		p1.wait()
		cmd = seal_base + " in=rand_sample.fq overwrite=t out=stdout.fq forbidn=f ambiguous=all stats=seal.stats rename=t trd=t ref={} hdist={} k={} t={}".format(args.ref, args.d, args.k, args.t)
	else:
		cmd = seal_base + " in={} overwrite=t out=stdout.fq forbidn=f ambiguous=all stats=seal.stats rename=t trd=t ref={} hdist={} k={} t={}".format(args.reads, args.ref, args.d, args.k, args.t)
		
	if args.m:
		cmd = cmd + " mm=t"
	else:
		cmd = cmd + " mm=f"

	sys.stderr.write("CMD: " + cmd + "\n")
	ferr = open("seal.log", "w")
	p1 = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=ferr)

	# seal is producing fastq output with the name line altered to include all of the
	# references hit by the read
	
	lnum = 0
	for szl in p1.stdout:
		szl = szl.strip()
		lnum += 1
		if lnum == 1:
			# name line
			# read name. parse out the targets
			aln = szl.split("\t")[1:]
			stats['aligned'] += 1

		elif lnum == 2:
			# read line, get the length of the read
			rlen = len(szl)
			
			# figure out number of kmers necessary to reach minimum percent aligned
			perc_aln = 0
			min_len = round(rlen*args.perc_aln)
			min_len_k = min_len - args.k + 1
			if min_len_k < 0:
				min_len_k = 0

			# handle it
			rnk = rlen - args.k + 1

			if len(aln) > 1:
				# potential multiple targets. parse out and sort by ratios
				targets = []
				for i in range(len(aln)):
					tmp = aln[i].split("=")
					kratio = float(tmp[1])/rnk
					if kratio >= args.min_ratio and int(tmp[1]) >= min_len_k and int(tmp[1]) >= args.min_kmers:
						targets.append([tmp[0], kratio])

				if len(targets) > 1:
					stats['assigned'] += 1

					# remaining targets after kratio filter
					if args.best:
						targets.sort(key=lambda x:x[1])
						targets.reverse()

						rbest = targets[0][1]
						i = 1
						while targets[i][1] == rbest:
							i += 1
							if i==len(targets):
								break
					else:
						i = len(targets)

					# if there were ties for best mapping then i is > 1
					bfactor = 1
					if i > 1:

						hfactor = 1.0/len(targets)

						if use_ambig:
							ghits = {}
							for hit in targets:
								ghits[ttog[hit[0]]] = 0

							if len(ghits.keys()) > 1:
								stats['gene_ambig'] += 1
								bfactor = 1.0/(float(len(ghits.keys()))**2)

						for hit in targets:
							hits[hit[0]][2] += hfactor
							hits[hit[0]][3] += hfactor*bfactor

						stats['ambig'] += 1

					else:
						# unique hit so assign it
						hits[targets[0][0]][1] += 1
						hits[targets[0][0]][2] += 1
						hits[targets[0][0]][3] += 1
						stats['unique'] += 1

				elif len(targets) == 1:
					# unique hit
					hits[targets[0][0]][1] += 1
					hits[targets[0][0]][2] += 1
					hits[targets[0][0]][3] += 1
					stats['assigned'] += 1
					stats['unique'] += 1
				else:
					# no hit, do nothing
					pass


			elif len(aln)==1:
				# only one target, assign it
				tmp = aln[0].split("=")
				hits[tmp[0]][1] += 1
				hits[tmp[0]][2] += 1
				hits[tmp[0]][3] += 1
				stats['assigned'] += 1
				stats['unique'] += 1

		elif lnum == 4:
			lnum = 0

	
	# finished parsing alignments
	p1.stdout.close()
	
	if stats['aligned'] > 0:

		# print stats
		sys.stderr.write("""
processed {} total reads (fragments) aligned; of these:
{} ({:0.2f}%) were assigned:
  unique:       {} ({:0.2f}%)
  ambiguous:    {} ({:0.2f}%)
  multi-gene:   {} ({:0.2f}%)

""".format(stats['aligned'], 
		stats['assigned'], float(stats['assigned'])/stats['aligned']*100,  
		stats['unique'], float(stats['unique'])/stats['assigned']*100,
		stats['ambig'], float(stats['ambig'])/stats['assigned']*100,
		stats['gene_ambig'], float(stats['gene_ambig'])/stats['assigned']*100))

	else:
		sys.stderr.write("\nNo alignments!\n")	

	# produce output
	print "\t".join(["transcript_id", "bundle", "length", "unique_hits", "hits", "adjusted_hits"])
	for tid in sorted(hits.keys()):
		bid = "-"
		if use_ambig:
			bid = ttog[tid]

		lout = [tid, bid] + hits[tid]
		print "\t".join(map(str, lout))
#		print "\t".join([
#					tid, bid, str(hits[tid][0]), 
#					str(hits[tid][1]), 
#					"{:0.4f}".format(hits[tid][2])])
	
	return 0





#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(
	description="""Naive method of expression quantification using seal (of bbmap) to
                   match reads via kmers to a reference""")

parser.add_argument("ref", type=str, 
	help="FASTA reference")
parser.add_argument("reads", type=str, 
	help="FASTQ or FASTA reads to align and quantify")

parser.add_argument("-r", dest="rand_sample", type=float, action="store", default=-1, 
	help="Random sample reads to this final number of reads (INT) or ratio of reads (float between 0 and 1) (default: all reads)")

# quantification options
parser.add_argument("-B", dest="B", action="store", type=str, default="",
	help="Transcript bundle file. A 2-column file with transcript ids and bundle ids. May be the result of a sequence clustering of transcripts into loci or genes.")

parser.add_argument("--min-ratio", dest="min_ratio", type=float, action="store", default=0, 
	help="Minimum ratio of kmers to accept as a valid assignment [0]")
parser.add_argument("--min-kmers", dest="min_kmers", type=int, action="store", default=1,
	help="Minimum number of kmers to accept as a match. [1]")
parser.add_argument("--perc-aln", dest="perc_aln", type=float, action="store", default=0, 
	help="Minimum percent of read (based on numer of kmers) for accepted hits [0]")
parser.add_argument("--best", dest="best", action="store_const", const=True, default=False, 
	help="Enable best mode which only assigns hits to the best match (the one with the most kmers). If there are multiple bests then the hit is divided between the targets.")

# seal options
parser.add_argument("-t", dest="t", action="store", type=int, default=1, 
	help="number of threads [1]")
parser.add_argument("-k", dest="k", action="store", type=int, default=15,
	help="kmer size [15]")
parser.add_argument("-d", dest="d", action="store", type=int, default=0,
	help="sets hdist parameter [0]")
parser.add_argument("-m", dest="m", action="store_const", const=True, default=False, 
	help="enables mm=t to allow a 'free' mismatch when matching kmers [off]")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

