#!/usr/bin/python
#
# juncs_to_psi.py
#
# Shawn Driscoll
# 20140618
#
# Parses output of juncs-and-boundaries.pl to produce psi/theta splicing 
# metrics. A graph is used to cluster junctions although this slows things down a bit.
#
# citation:
# Pervouchine DD1, Knowles DG, Guigo R. Intron-centric estimation of 
# alternative splicing from RNA-seq data. Bioinformatics. 2013 Jan 15;29(2):273-4
#

import sys, os
import igraph as ig
import argparse

def main(args):

	# set variables

	fname = args.juncs
	h5p = {}
	h3p = {}
	k1 = ""
	k2 = ""
	ivertex = 0
	v1 = 0
	v2 = 0
	juncs = []
	edges = []
	edge_counts = []
	v = []
	idx = 0

	# start a graph
	g = ig.Graph()

	# read junctions file

	sys.stderr.write("parsing {}\n".format(fname))
	fin = open(fname, "r")

	if args.s > 0:
		cnt = 0
		while cnt < args.s:
			szl = fin.readline()
			cnt += 1

	for szl in fin:

		arl = szl.strip().split("\t")

		# put all of the 5' and 3' ends into a hash with their overlap counts.
		# simultaneously count them and add their vertex ids. we'll also build a 
		# list of the junctions for later with the junction read depths.

		k1 = ":".join(arl[0:2])
		k2 = ":".join([arl[0], arl[2]])

		# add in ids to the donor/acceptor hashes
		if k1 not in h5p:
			h5p[k1] = [ivertex, 0, float(arl[9])]
			ivertex += 1
#			v.append(k1)

		if k2 not in h3p:
			h3p[k2] = [ivertex, 0, float(arl[10])]
			ivertex += 1
#			v.append(k2)

		juncs.append([k1, k2, float(arl[3])])

		# increment donor and acceptor depths
		h5p[k1][1] += float(arl[3])
		h3p[k2][1] += float(arl[3])

		# setup the edge for the graph. igraph uses integers for the vertices
		# so the h5p and h3p dicts are used to maintain the association between
		# vertex name and index.
		v1 = h5p[k1][0]
		v2 = h3p[k2][0]
		edges.append((v1, v2))
#		edge_counts.append(float(arl[3]))

	# data is loaded now we can play

	# build the graph
	g.add_vertices(ivertex)
	# add the edges
	g.add_edges(edges)
	# set the vertex names
#	g.vs["name"] = v
#	g.es["count"] = edge_counts
	# cluster junctions
	sys.stderr.write("clustering junctions and sorting by group\n")
	gi = g.clusters()

	v_to_g = gi.membership

	# add in the group ids to juncs and then sort it
	for i in range(len(juncs)):
		group = v_to_g[h5p[juncs[i][0]][0]]
		juncs[i].append(group)

	juncs.sort(key=lambda x: x[3])

	sys.stderr.write("writing results\n")

	for ji in juncs:
		v1 = ji[0]
		v2 = ji[1]
		count = float(ji[2])

		# get the cluster
		j = v_to_g[h5p[v1][0]]

		psi5 = count*1.0/h5p[v1][1]
		psi3 = count*1.0/h3p[v2][1]
		theta5 = h5p[v1][1]*1.0/(h5p[v1][1]+h5p[v1][2])
		theta3 = h3p[v2][1]*1.0/(h3p[v2][1]+h3p[v2][2])

		tmp = v1.split(":")
		tmp2 = v2.split(":")

		lout = [tmp[0], tmp[1], tmp2[1], "SC_{:08d}".format(j+1), count, h5p[v1][2], h3p[v2][2], psi5, psi3, theta5, theta3]
		print "\t".join(map(str, lout))		

	if False:

		sys.stderr.write("getting subgraphs\n")
		gc = gi.subgraphs()
		sys.stderr.write("generating results\n")

		for j in range(n):

			if j % 1000 == 0:
				sys.stderr.write("processed {} of {} clusters ({:0.2f}%)\n".format(j, n, j*100.0/n))

			gbar = gc[j]
			gbar_edges = gbar.get_edgelist()

			for i in range(len(gbar_edges)):
				edge = gbar_edges[i]
				v1 = gbar.vs[edge[0]]["name"]
				v2 = gbar.vs[edge[1]]["name"]

				if v1 not in h5p:
					# vertices are backwards because graph isn't directed
					tmp = v1
					v1 = v2
					v2 = tmp

				count = gbar.es[i]["count"]
				psi5 = count*1.0/h5p[v1][1]
				psi3 = count*1.0/h3p[v2][1]
				theta5 = h5p[v1][1]*1.0/(h5p[v1][1]+h5p[v1][2])
				theta3 = h3p[v2][1]*1.0/(h3p[v2][1]+h3p[v2][2])

				tmp = v1.split(":")
				tmp2 = v2.split(":")

				lout = [tmp[0], tmp[1], tmp2[1], "CLUST_{:08d}".format(j+1), count, h5p[v1][2], h3p[v2][2], psi5, psi3, theta5, theta3]
				print "\t".join(map(str, lout))

		sys.stderr.write("processed {} of {} clusters ({:0.2f}%)\n".format(n, n, n*100.0/n))


#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Calculate psi/theta splicing metrics from junctions output")
parser.add_argument("juncs", type=str, help="Juncs file with boundary counts (from juncs-and-boundaries.pl)")
parser.add_argument("-s", type=int, default=0, 
	help="Skip this many lines from top of file [0]")
args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")


