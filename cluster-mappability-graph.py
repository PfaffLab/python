#!/usr/bin/env python
#==============================================================================
# cluster-mappability-graph.py
#
# Shawn Driscoll
# 20141222
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Clusters the output of seal-mappability-graph or make-mapability-graph
#==============================================================================

import sys, argparse, math, re
import igraph as ig
import numpy as np
from scipy.cluster.vq import vq, kmeans, whiten

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
	g = ig.Graph()
	dverts = {}
	edges = []
	sims = []
	vidx = 0
	vnames = []

	# -- NOTE
	# sims and vnames aren't necessary since we only want to use the graph
	# to cluster features that have edge connections. 

	# open graph file and parse in the edges plus build vertices
	# dict
	fin = open(args.gin, "r")
	for szl in fin:
		aln = szl.strip().split("\t")
		if aln[0] not in dverts:
			dverts[aln[0]] = vidx
			vidx += 1
			vnames.append(aln[0])
		if aln[1] not in dverts:
			dverts[aln[1]] = vidx
			vidx += 1
			vnames.append(aln[1])

		# append edge tuple to list of edges
#		if float(aln[2]) >= args.t:
		edges.append((dverts[aln[0]], dverts[aln[1]]))
		sims.append(float(aln[2]))


	fin.close()

	# make a new set of edges that have sim scores passing the threshold
	# or in clusters determined via kmeans
	if args.k:
		sim_keeps = cluster_sims(sims)
	else:
		sim_keeps = list(sims)
		for i in range(len(sims)):
			if sims[i] >= args.t:
				sim_keeps[i] = 1
			else:
				sim_keeps[i] = 0

	edges_final = []
	for i in range(len(sims)):
		if sim_keeps[i] == 1:
			edges_final.append(edges[i])

	g.add_vertices(vidx)
	g.add_edges(edges_final)
	# add vertex names?
	g.vs["name"] = vnames
	# add similarity scores?
	g.es["sim"] = sims

	g_clusts = g.clusters()

	bundles = g_clusts.membership

	gout = []
	for i in range(len(bundles)):
		gout.append([vnames[i], bundles[i]])
	
	gout.sort(key=lambda x:x[1])

	for i in range(len(gout)):
		cid = "CLUST_{:08d}".format(gout[i][1])
		print "\t".join(map(str, [gout[i][0], cid]))

	return 0

def cluster_sims(v):
	# convert to numpy array
	v0 = np.array(v)
	# cluster
	cents, vss = kmeans(v0, 3)
	code, dists = vq(v0, cents)
	
	# generate cluster means
	gmeans = [0, 0, 0]
	gcounts = [0, 0, 0]
	for i in range(len(code)):
		gcounts[code[i]] += 1
		gmeans[code[i]] += v0[i]

	for i in range(3):
		gmeans[i] = gmeans[i]*1.0/gcounts[i]

	# figure out which is the minimum
#	mmin = 1
#	imin = 0
#	for i in range(3):
#		if gmeans[i] < mmin:
#			imin = i
#			mmin = gmeans[i]
#
#	# make a new cluster vector that has 1's where the middle and top
#	# mean groups were
#	gfinal = [0 for i in range(len(code))]
#	for i in range(len(code)):
#		if code[i] == imin:
#			gfinal[i] = 0
#		else:
#			gfinal[i] = 1
#
	# figure out which is the maximum
	mmax = 0
	imax = 0
	for i in range(3):
		if gmeans[i] > mmax:
			imax = i
			mmax = gmeans[i]

	# make a new cluster vector that has 1's where the top
	# mean group is
	gfinal = [0 for i in range(len(code))]
	for i in range(len(code)):
		if code[i] == imax:
			gfinal[i] = 1
		else:
			gfinal[i] = 0

	return(gfinal)

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
parser.add_argument('gin', type=str, help="File containing edges.")

parser.add_argument("-t", type=float, action="store", default=0, 
	help="Similarity threshold for edges [0]")
parser.add_argument("-k", action="store_const", const=True, default=False,
	help="Use clustering to determine edges to keep [off]")
parser.add_argument("-p", type=str, action="store", default="CLUST", 
	help="Prefix for cluster names [CLUST]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

