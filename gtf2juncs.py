#!/usr/bin/env python
#==============================================================================
# gtf2juncs.py
#
# Shawn Driscoll
# 20160819
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Print out a file of junctions from a GTF annotation. this will include 
# gene name, id, strand and transcript id info
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
	tid = ""
	gid = ""
	gname = ""
	ltid = ""
	laln = []
	jid = ""
	iidx = 0

	#
	# Start parsing the GTF
	#

	try:
		fin = open(args.gtf, "r")
	except:
		sys.stderr.write("Failed to open GTF ({})\n".format(args.gtf))
		return 1

	sys.stderr.write("Parsing {}\n".format(args.gtf))

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
				juncdb[jid] = dict(ref=aln[0], start=int(laln[4])-1, end=int(aln[3])-1, 
					strand=aln[6], gnames={}, gids={}, tids={})

			# add transcript id to the sets
			juncdb[jid]["tids"][tid] = 0

			# add gene id to the sets
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

		laln = list(aln)
		ltid = tid

	fin.close()

	# last line
	if tid == ltid:

		# junction is from laln[4] to aln[3]
		jid = "{}:{}-{}".format(aln[0], str(int(laln[4])+1), str(int(aln[3])-1))

		# is this junction in the database already?
		if jid not in juncdb:
			juncdb[jid] = dict(ref=aln[0], start=int(laln[4])-1, end=int(aln[3])-1, 
				strand=aln[6], gnames={}, gids={}, tids={})

		# add transcript id to the sets
		juncdb[jid]["tids"][tid] = 0

		# add gene id to the sets
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

	sys.stderr.write("done.\n")

	#
	# print bed type format
	#

	for jid in sorted(juncdb.keys()):

		# check length
		flen = juncdb[jid]['end']-juncdb[jid]['start']
		if flen < args.min_length or flen > args.max_length:
			# skip this one
			if flen < args.min_length:
				sys.stderr.write("Warning: skipping short intron ({}): {}\n".format(flen, jid))
			else:
				sys.stderr.write("Warning: skipping long intron ({}): {}\n".format(flen, jid))
			continue

		if args.g:
			iidx += 1
			iid = "INTID_{:08d}".format(iidx)

			lout = [
				juncdb[jid]['ref'], 
				"introns", "exon", 
				juncdb[jid]['start']+1, 
				juncdb[jid]['end'], 
				".", 
				juncdb[jid]['strand'], 
				".", 
				"gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; locus \"{}\"; oId \"{}\";".format(
					jid, iid, ",".join(juncdb[jid]['gnames']), iid, ",".join(juncdb[jid]['tids']))

			]
		else:
			lout = [
				juncdb[jid]['ref'], 
				juncdb[jid]['start'], 
				juncdb[jid]['end'], 
				jid, 
				juncdb[jid]['strand'],
				",".join(juncdb[jid]['gnames']), 
				",".join(juncdb[jid]['gids']),
				",".join(juncdb[jid]['tids'])]

		print "\t".join(map(str, lout))

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

#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="...")
parser.add_argument('gtf', type=str, help="GTF to extract junctions from")
parser.add_argument('-g', action="store_const", const=True, default=False, 
	help="Output in GTF format instead of BED.")
parser.add_argument('--min-length', default=10, type=int, action="store", 
	help="Minimum length of intron/junction feature to be exported [10]")
parser.add_argument('--max-length', default=200000, type=int, action="store", 
	help="Maximum length of intron/junction feature to be exported [200000]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

