#!/usr/bin/env python
#==============================================================================
# jucdb-quant-to-psi.py
#
# Shawn Driscoll
# 20150102
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys, argparse, math, re
import igraph as ig

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
	l3p = {}
	l5p = {}
	jid2p = {}
	dverts = {}
	lverts = []
	edges = []
	vidx = 0

	# check input file
	if not file_exists(args.infile):
		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.infile)
		return 1


	fin = open(args.infile, "r")
	for szl in fin:
		if re.search("unique_hits", szl):
			# make a dict of the column names
			aln = szl.strip().split("\t")
			lheader = {}
			for i in range(len(aln)):
				lheader[aln[i]] = i

			continue

		aln = szl.strip().split("\t")

		if re.search("[35]p$", aln[lheader["transcript_id"]]):
			# this is a specific donor/acceptor site and not a junction
			tid = aln[lheader['transcript_id']]
			if re.search("3p$", tid):
				# 3p site
				if tid not in l3p:
					l3p[id3p] = dict(juncs=[], junc_hits=[], total_hits=0, ovl_hits=0)

				l3p[id3p]["ovl_hits"] += int(aln[lheader['unique_hits']])

			if re.search("5p$", tid):
				# 5p site
				if tid not in l5p:
					l5p[id5p] = dict(juncs=[], junc_hits=[], total_hits=0, ovl_hits=0)

				l5p[id5p]["ovl_hits"] += int(aln[lheader['unique_hits']])

		else:

			# break up the junction id
			tmp = aln[lheader['transcript_id']].split(":")
			ref = tmp[0]
			tmp2 = tmp[1].split("-")
			pos3p = tmp2[1]
			pos5p = tmp2[0]

			id3p = "{}:{}:3p".format(ref, pos3p)
			id5p = "{}:{}:5p".format(ref, pos5p)

			# add vertices for graph of splice positions
			if id5p not in dverts:
				dverts[id5p] = vidx
				lverts.append(id5p)
				vidx += 1

			if id3p not in dverts:
				dverts[id3p] = vidx
				lverts.append(id3p)
				vidx += 1

			# add edge
			edges.append((dverts[id5p], dverts[id3p]))

			jid2p[aln[lheader['transcript_id']]] = [ref, pos5p, pos3p, aln[lheader['unique_hits']], id5p, id3p]

			if id3p not in l3p:
				l3p[id3p] = dict(juncs=[], junc_hits=[], total_hits=0, ovl_hits=0)

			if id5p not in l5p:
				l5p[id5p] = dict(juncs=[], junc_hits=[], total_hits=0, ovl_hits=0)

			l3p[id3p]["juncs"].append(aln[lheader['transcript_id']])
			l3p[id3p]['junc_hits'].append(aln[lheader['unique_hits']])
			l3p[id3p]['total_hits'] += int(aln[lheader['unique_hits']])

			l5p[id5p]["juncs"].append(aln[lheader['transcript_id']])
			l5p[id5p]['junc_hits'].append(aln[lheader['unique_hits']])
			l5p[id5p]['total_hits'] += int(aln[lheader['unique_hits']])

	fin.close()

	# cluster the nodes in the graph
	ngraph = ig.Graph()
	ngraph.add_vertices(range(len(lverts)))
	ngraph.add_edges(edges)
#	ngraph.vs['jname'] = lv2j

	gclusts = ngraph.clusters()
	vmem = gclusts.membership

#	# sort idx by the vertex clustering from the graph
#	idx = [x for (y, x) in sorted(zip(gclusts.membership, range(len(lverts))), key=lambda pair: pair[0])]

#	for i in idx:
#		print i, vmem[i], lverts[i]

	# make a dict that assocates junctions with cluster ids from the graph
	dtmp = {}

	for i in range(len(lverts)):
		if lverts[i] in l3p:
			for jid in l3p[lverts[i]]['juncs']:
				dtmp[jid] = vmem[i]
		elif lverts[i] in l5p:
			for jid in l5p[lverts[i]]['juncs']:
				dtmp[jid] = vmem[i]

	# change this to a list
	dtmp2 = []
	for tid in dtmp.keys():
		dtmp2.append((tid, dtmp[tid]))
	# sort this by the second value
	# jid_sorted = [x for (x, y) in sorted(dtmp2, key=lambda pair: pair[1])]

	print "\t".join(["chrom", "start", "end", "count", "locus", "donor_ovl", "donor_use", "acc_ovl", "acc_use", "psi5", "psi3", "theta5", "theta3"])

	# loop through the junctions and print everything out
	for jid in sorted(jid2p.keys()):

		juse = int(jid2p[jid][3])
		duse = int(l5p[jid2p[jid][4]]["total_hits"])
		dovl = int(l5p[jid2p[jid][4]]["ovl_hits"])
		ause = int(l3p[jid2p[jid][5]]["total_hits"])
		aovl = int(l3p[jid2p[jid][5]]["ovl_hits"])

		lout = [jid2p[jid][0], jid2p[jid][1], jid2p[jid][2], str(juse), "JLOC_{:08d}".format(dtmp[jid]), 
			str(dovl), str(duse), str(aovl), str(ause)]

		# psi5
		if duse > 0:
			lout.append("{:0.4f}".format(juse*1.0/duse))
		else:
			lout.append("NA")

		# psi3
		if ause > 0:
			lout.append("{:0.4f}".format(juse*1.0/ause))
		else:
			lout.append("NA")

		# theta5
		if (juse+dovl) > 0:
			lout.append("{:0.4f}".format(juse*1.0/(juse+dovl)))
		else:
			lout.append("NA")

		# theta3
		if (juse+aovl) > 0:
			lout.append("{:0.4f}".format(juse*1.0/(juse+aovl)))
		else:
			lout.append("NA")


		print "\t".join(map(str, lout))
		

	return 0


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
parser.add_argument('infile', type=str, help="Input file")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

