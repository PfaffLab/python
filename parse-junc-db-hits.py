#!/usr/bin/env python
#==============================================================================
# template.py
#
# Shawn Driscoll
# 20151202
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Pipeline to count mer hits to psi and theta splice site sequences. "alignments"
# come from Seal and input will be a fasta file containing the targets hit. 
#
# Updated 2016 to allow matching to a filtered MER index (juncdb-filter-matched-ref-mers.py)
#==============================================================================

import sys, argparse, math, re
import subprocess as sp

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

	# various stats variables
	num_alignments = 0
	num_alignments_counted = 0
	gids_hit = set()
	features_hit = set()
	features_counted = set()
	gids_counted = set()
	ambiguous_gid = 0

	if args.b:
		# if -b then we need -d too
		args.d = True

	# load GTF annotation
	sys.stderr.write("loading juncdb GTF annotation...\n")
	annot = parse_gtf(args.gtf)
	sys.stderr.write("found {} features\n".format(len(annot.keys())))

	# donor and acceptor tables. each will contain total usage and 
	# total read-through
	tdonor = {}
	tacc = {}
	sys.stderr.write("building donor and acceptor tables from annotation...\n")
	# build out the donor and acceptor tables
	for tid in annot.keys():
		# check type
#		if annot[tid]['type'] == "theta":
#
#			if re.search("5p$", tid):
#				# this is a donor
#				tdonor[tid] = {'use': 0, 'pass': 0}
#			else:
#				tacc[tid] = {'use': 0, 'pass': 0}

		if annot[tid]['type'] == "psi":
			# this is a junction which has a donor and an acceptor...perfect!
			rres = re.search("^([^\:]+)\:([0-9]+)\-([0-9]+)", tid)
			if rres:
				lside = ":".join([rres.group(1), rres.group(2)])
				rside = ":".join([rres.group(1), rres.group(3)])

				if annot[tid]['strand'] == "+":
					tdonor[lside] = {'use': 0, 'pass': 0}
					tacc[rside] = {'use': 0, 'pass': 0}
				else:
					tdonor[rside] = {'use': 0, 'pass': 0}
					tacc[lside] = {'use': 0, 'pass': 0}

			else:
				sys.stderr.write("Problem parsing junction id {}\n".format(tid))


	# now open up the fasta and get serious
	rres = re.search("\\.gz$", args.fasta)
	if rres:
		sys.stderr.write("FASTA file is gzipped. I can deal with that...\n")
		p1 = sp.Popen("gunzip -c {}".format(args.fasta).split(), stdout=sp.PIPE)
	else:
		p1 = sp.Popen("cat {}".format(args.fasta).split(), stdout=sp.PIPE)

	sys.stderr.write("assigning hits...\n")
	for szl in p1.stdout:
		rres = re.search("^>", szl)
		if rres:
			num_alignments += 1
			
			# name row, parse this for hits
			szl = re.sub('^>', '', szl)
			aln = szl.strip().split("\t")

			# track gene ids hit
			tgid = {}

			# expand the targets
			targets = []
			targets_raw_hits = []

			# update 12 Feb 2016
			# now we are matching reads to already kmer'd mers so each read may 
			# return hits to the same feature but dozens of times so the first pass 
			# has to be to bin them and sum their counts. this next loop will basically
			# just reformat the new input to what the old input looked like so the
			# remaining code can just carry on like nothing happened
			#

			tbin = {}
			for i in range(1, len(aln)):
				atmp = aln[i].split("=")
				atmp2 = atmp[1].split("|")

				for tid in atmp2:
					if tid not in tbin:
						tbin[tid] = 0

					tbin[tid] += int(atmp[2])

			aln0 = list(aln)
			aln = [aln0[0]]
			for tid in tbin.keys():
				aln.append("=".join([tid, str(tbin[tid])]))
			

			for i in range(1, len(aln)):

				#
				# the name consists of the 'tid' key for the annot dict as well as a hit count
				# that follows an '=' sign. so we can split on that. if we are counting hits
				# from a kmer database (i.e args.m==True) then the hit count is the feature 
				# name in the FASTA line. If we are couting hits from raw FASTQ then the
				# kmer hit count is in each of these target/count pairs. If the count unit
				# is reads insted of kmers then each count is one unless modified below.
				#

				# split the hit
				atmp = aln[i].split("=")
				if args.m:
					# name of the fasta row is the count because it was a kmer database
					atmp[1] = int(aln[0])
				else:
					# the number of kmer matches is in atmp[1]
					targets_raw_hits.append(int(atmp[1]))

					if int(atmp[1]) >= args.n:
						if args.k:
							atmp[1] = int(atmp[1])
						else:
							atmp[1] = 1
					else:
						atmp[1] = 0

				# update set of all features hit by the entire file
				features_hit.update([atmp[0]])

				# check gene id of the hit and append it
				if annot[atmp[0]]['gid'] not in tgid:
					tgid[annot[atmp[0]]['gid']] = []

				# append the index of this target in the targets list
				tgid[annot[atmp[0]]['gid']].append(i-1)

				# update global counter of gene features hit
				gids_hit.update([annot[atmp[0]]['gid']])

				# check the type of this feature - psi or theta
				ttype = annot[atmp[0]]['type']
			 	
				if ttype=="psi":
					# psi hit, update 
					lside, rside = psi_name_to_da(atmp[0])
					atmp.append("psi")

					# look at strand and append as donor, acceptor
					if annot[atmp[0]]['strand'] == "+":
						atmp.append(lside)
						atmp.append(rside)
					else:
						atmp.append(rside)
						atmp.append(lside)

				else:
					atmp.append("theta")
					lside = theta_name_to_da(atmp[0])
					atmp.append(lside)

				targets.append(list(atmp))

			if args.d and args.b and (not args.m):
				# in this mode if there are multiple gene ids hit by this read then
				# we have to figure out if there is a 'best' match and discard those
				# with less than best matches
				
				if len(targets) > 1 and len(tgid.keys()) > 1:
					
					tgid_hits = {}
					gid_raw_hits = []
					tgid_names = tgid.keys()

					for gid in tgid_names:
						tgid_hits[gid] = 0
						for k in range(len(tgid[gid])):
							tgid_hits[gid] += targets_raw_hits[k]

					for gid in tgid_names:
						gid_raw_hits.append(tgid_hits[gid])

					# sort hit indices in descending order
					o = list(order(gid_raw_hits))
					o.reverse()

					targets_hat = []
					tgid_hat = {}
					for i in range(len(tgid_names)):
						if i==0:
							# append this one's targets
							for j in tgid[tgid_names[o[0]]]:
								targets_hat.append(list(targets[j]))

							tgid_hat[tgid_names[o[0]]] = 0

						else:
							if gid_raw_hits[o[i]] == gid_raw_hits[o[0]]:
								# keep this one too
								for j in tgid[tgid_names[o[i]]]:
									targets_hat.append(list(targets[j]))

								tgid_hat[tgid_names[o[i]]] = 0

					targets = targets_hat
					tgid = tgid_hat
			
			elif not args.d:
				# hits to multiple gene ids are not allowed in this mode
				if len(tgid.keys()) > 1:
					tgid = {}
					targets = []

			if len(targets) > 0:

				#
				# if this feature hits both normal junctions and either of the
				# two retention types then it is ambiguous and we don't want to 
				# count it at all especially since I supposedly removed all of 
				# those kind of mers when making the reference
				#

				hit_junc = False
				hit_ret3 = False
				hit_ret5 = False
				for i in range(len(targets)):
					rres = re.search("5p$", targets[i][0])
					if re.search("5p$", targets[i][0]):
						hit_ret5 = True

					if re.search("3p$", targets[i][0]):
						hit_ret3 = True

					if re.search("[0-9]+\-[0-9]+$", targets[i][0]):
						hit_junc = True

				if ((hit_ret3 or hit_ret5) and hit_junc) or (hit_ret3 and hit_ret5):
					# throw this out
					targets = []
					tgid = {}


			if len(tgid.keys()) > 0:

				num_alignments_counted += 1
				for gid in tgid.keys():
					gids_counted.update([gid])

				if len(tgid.keys()) > 1:
					ambiguous_gid += 1

				# down-weighting for multiple gene ids is 1/(number of gids)^2
				tgid_weight = float(1)/(len(tgid.keys())**2)

				# loop through targets and assign hits
				for i in range(len(targets)):

					features_counted.update([targets[i][0]])

					hits = (targets[i][1]*tgid_weight)
					# first assign the hit to the annot table
					annot[targets[i][0]]['hits'] += hits

					if targets[i][2] == "psi":
						
						# junction hit, add total to donor and acceptor as well
						if targets[i][3] in tdonor:
							tdonor[targets[i][3]]['use'] += hits
						if targets[i][4] in tacc:
							tacc[targets[i][4]]['use'] += hits

					else:
						if targets[i][3] in tdonor:
							tdonor[targets[i][3]]['pass'] += hits
						elif targets[i][3] in tacc:
							tacc[targets[i][3]]['pass'] += hits
						else:
							sys.stderr.write("[warning] Couldn't find overlap site in the tables ({})\n".format(targets[i][3]))

			# print status to keep the user involved...
			if (num_alignments % 1000000) == 0:
				sys.stderr.write("parsed/assigned {}/{} reads\n".format(num_alignments, num_alignments_counted))

	# print out statistics...
	sys.stderr.write("\n")
	sys.stderr.write("Total reads/kmers parsed/counted:  {}/{}\n".format(num_alignments, num_alignments_counted))
	sys.stderr.write("Total genes hit/counted:           {}/{}\n".format(len(gids_hit), len(gids_counted)))
	sys.stderr.write("Total features hit/counted:        {}/{}\n".format(len(features_hit), len(features_counted)))
	sys.stderr.write("Total with ambiguous gene:         {}\n".format(ambiguous_gid))
	sys.stderr.write("\n")
	sys.stderr.write("Printing results to stdout...\n")

	#
	# hits are all assigned. now we can produce the output
	#

	# print header
	print "\t".join(["chrom", "start", "end", "count", "id", "strand", "donor_ovl", "donor_use", "acc_ovl", "acc_use", 
		"psi5", "psi3", "theta5", "theta3", "gene_name", "gene_id", "transcript_id"])

	for tid in sorted(annot.keys()):

		lout = []
		# only need to output at the level of junctions
		if annot[tid]['type'] == "psi":
			rres = re.search("^([^\:]+)\:([0-9]+)\-([0-9]+)", tid)
			if rres:
				lout = [rres.group(1), rres.group(2), rres.group(3), annot[tid]['hits'], tid]
				lside, rside = psi_name_to_da(tid)
				# get the counts for this thing
				# lout.append(annot[tid]['hits'])
				lout.append(annot[tid]['strand'])

				# get donor and acceptor use and overlap counts
				if annot[tid]['strand']=="+":
					did = lside
					aid = rside
				else:
					did = rside
					aid = lside

				dvals = [tdonor[did]['pass'], tdonor[did]['use']]
				avals = [tacc[aid]['pass'], tacc[aid]['use']]

				#
				# calculate metrics
				#

				# psi5, psi3, theta5, theta3
				metrics = ["NA", "NA", "NA", "NA"]

				if dvals[1] > 0:
					metrics[0] = float(annot[tid]['hits'])/dvals[1]
				if dvals[0] > 0 or dvals[1] > 0:
					metrics[2] = float(dvals[1])/(dvals[0]+dvals[1])

				if avals[1] > 0:
					metrics[1] = float(annot[tid]['hits'])/avals[1]
				if avals[0] > 0 or avals[1] > 0:
					metrics[3] = float(avals[1])/(avals[0]+avals[1])

				#
				# append counts and metrics
				#
				
				lout += dvals
				lout += avals
				lout += metrics

				#
				# append annotation
				#
				lout.append(annot[tid]['gname'])
				lout.append(annot[tid]['gid'])
				lout.append(",".join(annot[tid]['tid']))

				print "\t".join(map(str, lout))

	return 0

def psi_name_to_da(sz):

	lside = ""
	rside = ""

	# parse it apart with regex
	rres = re.search("^([^\:]+)\:([0-9]+)\-([0-9]+)", sz)
	if rres:
		lside = ":".join([rres.group(1), rres.group(2)])
		rside = ":".join([rres.group(1), rres.group(3)])
	else:
		sys.stderr.write("[warning] unable to parse psi name: {}\n".format(sz))

	return (lside, rside)

def theta_name_to_da(sz):

	atmp = sz.split(":")
	return ":".join(atmp[0:2])


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

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

		# parse attributes
		attr = parse_gtf_attr(aln[8])

		if attr['transcript_id'] not in gtfdb:
			gtfdb[attr['transcript_id']] = {'gid': attr['gene_id'], 'gname': attr['gene_name'], 
			'tid': delim_list_to_set(attr['oId'], ","), 'hits': 0, 'type': aln[1], 'strand': aln[6]}

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

def delim_list_to_set(sz, delim):
	aa = sz.split(delim)
	da = {}

	for n in aa:
		da[n] = 0

	return da

# like the R order function which returns a sorting of the indices
def order(v):
	return np.argsort(v)

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Parse the output of seal for counts of kmers to targets.")
parser.add_argument('fasta', type=str, help="Input file")
parser.add_argument('gtf', type=str, help="Junction database GTF (from make-junc-db.py)")
parser.add_argument('-m', action="store_const", const=True, default=False, 
	help="Use if what was counted in Seal was a kmer database (with counts as names) and not raw FASTQ reads")
parser.add_argument("-k", action="store_const", const=True, default=False,
	help="Use the kmer as the count unit instead of reads [off]")
parser.add_argument("-n", action="store", default=1, type=int, 
	help="Minimum kmer-count at a feature to count as a hit (not used with -m) [1]")
parser.add_argument("-b", action="store_const", const=True, default=False, 
	help="Unless using -m this option evaluates gene ids hit by a read and only counts hits to the one with the best match (most kmers) [off]")
parser.add_argument("-d", action="store_const", const=True, default=False, 
	help="Downweight multi-gene matched hits (default is to discard them) [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

