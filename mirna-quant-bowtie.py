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
import igraph as ig

def main(args):

	hits = {}
	hfactor = 0
	stats = dict(reads=0, aligned=0, assigned=0, unique=0, ambig=0, gene_ambig=0, lmm=0, proc=0)
	ghits = {} # used for gene/locus hit ambiguity counting
	ttog = {} # for translating transcript id to gene/locus id for ambiguity
	use_ambig = False
	dverts = {}
	edges = []
	sims = []
	vidx = 0
	vnames = []
	targets = []

	ambig_bundles = {}

	# bowtie2 command
	bowtie = "bowtie2 -p {} -k 100 --gbar 1000 --local -x {} -U {}"

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
	# run alignments
	# --
	
	cmd = bowtie.format(args.p, args.ref, args.reads)
	sys.stderr.write("CMD: " + cmd + "\n")
	ferr = open("bowtie2.log", "w")
	p1 = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=ferr)

	# --
	# parse alignments
	# --

	sys.stderr.write("parsing alignments from bowtie2...\n")

	cname = ""
	lname = ""
	vidx = 0

	for szl in p1.stdout:
		szl = szl.strip()

		if szl[0] == "@":
			aln = szl.split("\t")
			if aln[0] == "@SQ":
				tid = aln[1].split(":")[1]
				tlen = aln[2].split(":")[1]
				hits[tid] = [int(tlen), 0, 0]
				# add vertext
				dverts[tid] = vidx
				vidx += 1

		else:
			# alignent
			stats['proc'] += 1

			if stats['proc'] % 1e6 == 0:
				sys.stderr.write("reads: {}; aligned: {}; assigned: {}\n".format(stats['reads'], stats['aligned'], stats['assigned']))

			aln = szl.split("\t")
			cname = aln[0]

			if int(aln[1]) & 0x4:
				stats['reads'] += 1
				lname = cname
				continue

			# aligned
			if cname != lname:
				stats['reads'] += 1
				stats['aligned'] += 1

				if len(targets) > 0:
					# deal with alignments
					stats['assigned'] += 1

					# sort targets if necessary and keep only those that are "best"
					if len(targets) > 1:
						targets.sort(key=lambda x:x[1])
						targets.reverse()
						tmp = [targets[0]]
						i = 1
						while targets[i][1] == targets[0][1]:
							tmp.append(targets[i])
							i += 1
							if i == len(targets):
								break
						targets = list(tmp)

					if len(targets) == 1:
						# one target, done deal
						stats['unique'] += 1
						hits[targets[0][0]][1] += 1
						hits[targets[0][0]][2] += 1

					else:
						# ambiguous with multiple "best" alignments

						lbid = []
						for j in range(len(targets)):
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


				targets = []

			# append alignment to targets if it passes requirement
			rres = re.search("NM\:i\:([0-9]+)", szl)
			aln_nm = int(rres.group(1))
			if aln_nm <= 1+args.u:
				targets.append([aln[2], aln_nm])

			lname = cname


	sys.stderr.write("reads: {}; aligned: {}; assigned: {}\n".format(stats['reads'], stats['aligned'], stats['assigned']))

	# handle final read if any
	if len(targets) > 0:
		# deal with alignments
		stats['assigned'] += 1

		# sort targets if necessary and keep only those that are "best"
		if len(targets) > 1:
			targets.sort(key=lambda x:x[1])
			targets.reverse()
			if targets[0][1] <= 1:
				tmp = [targets[0]]
				i = 1
				while targets[i][1] <= (targets[0][1]+args.u):
					tmp.append(targets[i])
					i += 1
					if i == len(targets):
						break
				targets = list(tmp)
			else:
				# best alignment has more than 1 mismatch, drop the assignment
				targets = []

		if len(targets) == 1:
			# one target, done deal
			if targets[0][1] <= 1:
				stats['unique'] += 1
				hits[targets[0][0]][1] += 1
				hits[targets[0][0]][2] += 1
			else:
				stats['assigned'] -= 1

		elif len(targets) > 1:
			# ambiguous with multiple "best" alignments

			lbid = []
			for j in range(len(targets)):
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
			# no targets, must have had too many mismatches in best alignment
			stats['assigned'] -= 1


	# finished parsing alignments
	p1.stdout.close()

	if stats['assigned'] == 0:
		sys.stderr.write("No reads were assigned!\n")
		return(1)

	# --
	# now that the initial parsing and counting is done we have to deal with the 
	# ambiguous bundles. this is where things get a little interesting
	# --

	sys.stderr.write("assigning ambiguous hits\n")

	# for each bundle, get the mappability ratios from mapp and the # of 
	# unique hits for each target in the bundle. divide the unique hits 
	# by the ratios then divide those by their total sum. then the total hits
	# in the bundle is multiplied by those new ratios and assigned to the targets
	# summing in their unique hit counts
	for bid in ambig_bundles.keys():
		tid = bid.split("|")
		weights = [0 for i in range(len(tid))]
		wsum = 0
		for i in range(len(tid)):
			weights[i] = hits[tid[i]][1]*1.0 + 1
			wsum += weights[i]

		if wsum > 0:
			# loop back through and scale weights but sum so they sum to 1
			for i in range(len(tid)):
				weights[i] /= float(wsum)

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
processed {} total reads; of these {} ({:0.1f}%) were aligned and:
  {} ({:0.2f}%) were assigned:
    unique:       {} ({:0.2f}%)
    ambiguous:    {} ({:0.2f}%)
    multi-gene:   {} ({:0.2f}%)

""".format(stats['reads'], stats['aligned'], float(stats['aligned'])/stats['reads']*100, 
		stats['assigned'], float(stats['assigned'])/stats['aligned']*100,  
		stats['unique'], float(stats['unique'])/stats['assigned']*100,
		stats['ambig'], float(stats['ambig'])/stats['assigned']*100,
		stats['gene_ambig'], float(stats['gene_ambig'])/stats['assigned']*100))
	
	# produce output
	lout = []

	# get sum total hits
	total_hits = 0
	for tid in hits.keys():
		total_hits += hits[tid][2]

	total_rpkm = 0
	for tid in hits.keys():
		# calculate rpkm for features. keep a total sum and then we'll print out TPM
		# at the end as well
		rpkm = hits[tid][2]*1e9/(total_hits+hits[tid][0])
		total_rpkm += rpkm
		hits[tid].append(rpkm)

	print "\t".join(["transcript_id", "bundle", "length", "unique_hits", "hits", "rpkm", "tpm"])
	for tid in sorted(hits.keys()):
		bid = "-"
		if use_ambig:
			bid = "BID_{:08d}".format(ttog[tid])
		else:
			# use the graph bundles
			bid = "BID_{:08d}".format(bundles[dverts[tid]])
		print "\t".join([
					tid, bid,
					str(hits[tid][0]), 
					"{:0.4f}".format(hits[tid][1]),
					"{:0.4f}".format(hits[tid][2]),
					"{:0.4f}".format(hits[tid][3]),
					"{:0.4f}".format(hits[tid][3]*1e6/total_rpkm)])
	
	return 0




#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(
	description="aligns trimmed reads with bowtie2 and counts hits to features")

parser.add_argument("ref", type=str, 
	help="Bowtie2 index")
parser.add_argument("reads", type=str, 
	help="FASTQ")

# quantification options
parser.add_argument("-B", dest="B", action="store", type=str, default="",
	help="Transcript bundle file. A 2-column file with transcript ids and bundle ids. May be the result of a sequence clustering of transcripts into loci or genes.")

# seal options
parser.add_argument("-p", dest="p", action="store", type=int, default=1, 
	help="number of threads for bowtie and sorting [1]")
parser.add_argument("-u", dest="u", action="store", type=int, default=0,
	help="Number of mismatches in addition to the number in the best alignment to consider for assignment [0]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

