#!/usr/bin/python
#
# bwa-quant.py
#
# porting the perl version over to python...just because
#
#

import os,sys,re
import argparse
#from hashlib import md5
import subprocess as sp
import igraph as ig

def main(args):

	hits = {}
	hfactor = 0
	stats = dict(aligned=0, assigned=0, unique=0, ambig=0, gene_ambig=0)
	ghits = {} # used for gene/locus hit ambiguity counting
	ttog = {} # for translating transcript id to gene/locus id for ambiguity
	use_ambig = False
	dverts = {}
	edges = []
	sims = []
	vidx = 0
	vnames = []

	ambig_bundles = {}

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
		fin.close()

	# --
	# measure mappability of the reference
	# --
	sys.stderr.write("evaluating mappability of targets at k = {}\n".format(args.k))
	mapp = measure_mappability(args.ref, args.k, args.t)
	sys.stderr.write("finished\n".format(args.k))

	fout = open("mappability.txt", "w")
	for tid in sorted(mapp.keys()):
		fout.write("\t".join([tid, "{:0.4f}".format(mapp[tid])]) + "\n")
	fout.close()


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
		hits[aln[0]] = [int(aln[1]), 0, 0]
		dverts[aln[0]] = vidx
		vidx += 1

	fin.close()
	sys.stderr.write("finished\n".format(args.k))

	# -- 
	# run seal to start kmer matching process to the reference
	# --
	
#	sys.stderr.write("finished\n".format(args.k))

	# run alignment(s)
	cmd = seal_base + " in={} out=stdout.fq rename=t trd=t ref={} hdist={} k={} t={}".format(args.reads, args.ref, args.d, args.k, args.t)
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
	perc_aln = 0

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
				# multiple targets. parse out and sort by ratios
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
					if i > 1:
						# assign ambig hits to ambig bundles
						lbid = []
						for j in range(i):
							lbid.append(targets[j][0])
						
						lbid.sort(key=lambda x:x)
						bid = "|".join(lbid)
						if bid not in ambig_bundles:
							ambig_bundles[bid] = 0

						# add edges to edge list for graph
						for i in range(len(lbid)-1):
							for j in range(i, len(lbid)):
								edges.append((dverts[lbid[i]], dverts[lbid[j]]))

						ambig_bundles[bid] += 1	
						stats['ambig'] += 1

					else:
						# unique hit so assign it
						hits[targets[0][0]][1] += 1
						hits[targets[0][0]][2] += 1
						stats['unique'] += 1

				elif len(targets) == 1:
					# unique hit
					hits[targets[0][0]][1] += 1
					hits[targets[0][0]][2] += 1
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
				stats['assigned'] += 1
				stats['unique'] += 1

		elif lnum == 4:
			lnum = 0


	# finished parsing alignments
	p1.stdout.close()

	# --
	# now that the initial parsing and counting is done we have to deal with the 
	# ambiguous bundles. this is where things get a little interesting
	# --

	sys.stderr.write("assigning ambiguous hits\n".format(args.k))

	# for each bundle, get the mappability ratios from mapp and the # of 
	# unique hits for each target in the bundle. divide the unique hits 
	# by the ratios then divide those by their total sum. then the total hits
	# in the bundle is multiplied by those new ratios and assigned to the targets
	# summing in their unique hit counts
	for bid in sorted(ambig_bundles.keys()):
		tid = bid.split("|")
		weights = [0 for i in range(len(tid))]
		wsum = 0
		for i in range(len(tid)):
			if mapp[tid[i]] > 0:
				weights[i] = hits[tid[i]][1]*1.0/mapp[tid[i]];
			wsum += weights[i]

		if wsum > 0:
			# loop back through and scale weights but sum so they sum to 1
			for i in range(len(tid)):
				weights[i] /= wsum
		else:
			for i in range(len(tid)):
				weights[i] = 1.0/len(tid)

		# now assign the hits
		for i in range(len(tid)):
			hits[tid[i]][2] += ambig_bundles[bid]*weights[i]


	# cluster graph
	g = ig.Graph()
	g.add_vertices(len(dverts.keys()))
	g.add_edges(edges)
	g_clust = g.clusters()
	bundles = g_clust.membership

	# print stats
	sys.stderr.write("""
aligned {} total reads (fragments); of these:
{} ({:0.2f}%) were assigned:
  unique:       {} ({:0.2f}%)
  ambiguous:    {} ({:0.2f}%)
  multi-gene:   {} ({:0.2f}%)

""".format(stats['aligned'], 
		stats['assigned'], float(stats['assigned'])/stats['aligned']*100,  
		stats['unique'], float(stats['unique'])/stats['assigned']*100,
		stats['ambig'], float(stats['ambig'])/stats['assigned']*100,
		stats['gene_ambig'], float(stats['gene_ambig'])/stats['assigned']*100))
	

	# produce output
	print "\t".join(["transcript_id", "bundle", "length", "unique_hits", "hits"])
	for tid in sorted(hits.keys()):
		bid = "-"
		if use_ambig:
			bid = ttog[tid]
		else:
			# use the graph bundles
			bid = str(bundles[dverts[tid]])
		print "\t".join([
					tid, bid,
					str(hits[tid][0]), 
					"{:0.4f}".format(hits[tid][1]),
					"{:0.4f}".format(hits[tid][2])])
	
	return 0

# --
# measure_mappability
# run seal on the reference vs itself to measure mappability at the 
# selected kmer size
def measure_mappability(ref, k, t):
	seal_cmd = "seal.sh in={} ref={} k={} rename=t trd=t mm=f out=stdout.fa t={}"
	ferr = open("/dev/null", "w")
	p1 = sp.Popen(seal_cmd.format(ref, ref, k, t).split(), stdout=sp.PIPE, stderr=ferr)
	mapp = {}

	# parse the result
	for szl in p1.stdout:
		if szl[0] == ">":
			aln = szl.strip().split("\t")
			aln[0] = re.sub("^>", "", aln[0])
			hits = []
			refk = 0
			for i in range(1, len(aln)):
				tmp = aln[i].split("=")
				if tmp[0] == aln[0]:
					refk = int(tmp[1])
				else:
					hits.append([tmp[0], int(tmp[1])])

			# default ratio is 1.0 (100% mappable) if there were no other hits
			ratio = 1.0
			if len(hits) > 0:
				# there were other hits so the highest one is what we want
				hits.sort(key=lambda x:x[1])
				hits.reverse()
				# best one is the first one
				ratio = 1 - hits[0][1]*1.0/refk

			mapp[aln[0]] = ratio

	ferr.close()
	p1.stdout.close()

	# all done
	return(mapp)



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
parser.add_argument("-k", dest="k", action="store", type=int, default=13,
	help="kmer size [13]")
parser.add_argument("-d", dest="d", action="store", type=int, default=1,
	help="sets hdist parameter [1]")
parser.add_argument("-m", dest="m", action="store_const", const=True, default=False, 
	help="enables mm=t to allow a 'free' mismatch when matching kmers [off]")
parser.add_argument("--mp", dest="mapp", action="store", type=str, default="", 
	help="mappability vector. if not provided then the code makes one")


args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

