#!/usr/bin/env python
#==============================================================================
# bam2juncsA.py
#
# Shawn Driscoll
# 20160412
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# This script starts by finding junction alignments and building a junction 
# database with counts.
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
import hashlib
import subprocess as sp
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")
HBIN = 16000

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	djuncs = {}
	
	dp = {}

	id5p = ""
	id3p = ""

	jid = ""
	rcount = 0
	jcount = 0
	rjcount = 0
	lcount = 0

	# -------------------------------------------------------------------------
	# step 1: parse BAM and find all of the junctions
	# -------------------------------------------------------------------------

	# open bam

	if re.search("\.bam$", args.bam):
		p1 = sp.Popen("samtools view -F 0x4 {}".format(args.bam).split(), stdout=sp.PIPE)
	elif re.search("\.sam$", args.bam):
		p1 = sp.Popen("samtools view -SF 0x4 {}".format(args.bam).split(), stdout=sp.PIPE)

	for szl in p1.stdout:

		lcount += 1
		aln = parse_bam_alignment(szl.strip())
		if aln_is_aln(aln):
			rcount += 1

		# check cigar
		if re.search("[0-9]+N", aln["cigar"]):
			rjcount += 1
			lcig = explode_cigar(aln["cigar"])
			# step through to find the positions of the junction
			left0 = int(aln["pos"])
			left = left0

			for i in range(len(lcig)):
				if lcig[i][1]=="M" or lcig[i][1]=="D":
					left += int(lcig[i][0])
				elif lcig[i][1]=="N":
					right = left+int(lcig[i][0])-1
					jid = "{}:{}-{}".format(aln['rname'], left, right)
					id5p = "{}:{}:5".format(aln['rname'], left)
					id3p = "{}:{}:3".format(aln['rname'], right)

					if jid not in djuncs:
						jcount += 1
						djuncs[jid] = [aln['rname'], left, right, hashlib.md5(jid).hexdigest(), 0, 
							id5p, id3p]

					if args.psi_theta:
						if id5p not in dp:
							# two counts - total usage and total overlap
							dp[id5p] = [aln['rname'], left, 0, 0]
						if id3p not in dp:
							dp[id3p] = [aln['rname'], right, 0, 0]

						# increment counts for this 5' and 3' end
						dp[id5p][2] += 1
						dp[id3p][2] += 1

					# increment count for this junction
					djuncs[jid][4] += 1
					# move left point to the next cigar op
					left += int(lcig[i][0])

		if lcount % 1000000 == 0:
			sys.stderr.write("L:{}; A:{}({:0.2f}); JA:{}({:0.2f}); J:{}".format(lcount, rcount, 
				rcount*1.0/lcount, rjcount, rjcount*1.0/rcount, jcount) + "\n")

	sys.stderr.write("L:{}; A:{}({:0.2f}); JA:{}({:0.2f}); J:{}".format(lcount, rcount, 
		rcount*1.0/lcount, rjcount, rjcount*1.0/rcount, jcount) + "\n")
	
	if args.psi_theta:
		# ---------------------------------------------------------------------
		# step 2: read the BAM again but this time check for the overlaps at 
		# 5' and 3' ends of junctions
		# ---------------------------------------------------------------------

		# build the lookup hash for 5' and 3'
		sys.stderr.write("Building 5' and 3' lookup tables\n")

		# build list of the "primes"
		lp = []
		for idp in dp.keys():
			lp.append(idp.split(":"))

		hp = build_point_lookup_hash(lp, HBIN)

		if re.search("\.bam$", args.bam):
			p1 = sp.Popen("samtools view -F 0x4 {}".format(args.bam).split(), stdout=sp.PIPE)
		elif re.search("\.sam$", args.bam):
			p1 = sp.Popen("samtools view -SF 0x4 {}".format(args.bam).split(), stdout=sp.PIPE)

		sys.stderr.write("Checking for alignments overlapping 5' and 3' junction ends\n")

		lcount = 0
		rcount = 0
		rjcount = 0

		for szl in p1.stdout:
			lcount += 1
			aln = parse_bam_alignment(szl.strip())
			if not aln_is_aln(aln):
				continue

			rcount += 1
			mm = aln_get_mapped_regions(aln)
			for i in range(len(mm)):
				# check each one for hits to something
				hit = False
				tmp = find_point_hash_hits(mm[i], hp, HBIN)
				if len(tmp) > 0:
					hit = True
					# assign hits
					for pid in tmp:
						dp[pid][3] += 1

				if hit:
					rjcount += 1

			if lcount % 1000000 == 0:
				sys.stderr.write("L:{}; A:{}({:0.2f}); OA:{}({:0.2f})".format(lcount, rcount, 
					rcount*1.0/lcount, rjcount, rjcount*1.0/rcount) + "\n")

		sys.stderr.write("L:{}; A:{}({:0.2f}); OA:{}({:0.2f})".format(lcount, rcount, 
			rcount*1.0/lcount, rjcount, rjcount*1.0/rcount) + "\n")



	if args.psi_theta:

		for jid in sorted(djuncs.keys()):
			lout = djuncs[jid][0:5]
			id5p = djuncs[jid][5]
			id3p = djuncs[jid][6]

			# total 5p use, total 3p use, total 5p ovl, total 3p ovl, psi5, psi3, theta5, theta3
			psi5 = 0
			psi3 = 0
			theta5 = 0
			theta3 = 0
			u5p = dp[id5p][2]
			o5p = dp[id5p][3]
			u3p = dp[id3p][2]
			o3p = dp[id3p][3]

			if u5p > 0:
				psi5 = lout[4]*1.0/u5p
				theta5 = u5p*1.0/(u5p+o5p)
			else:
				psi5 = "NA"
				if o5p > 0:
					theta5 = 0
				else:
					theta5 = "NA"

			if u3p > 0:
				psi3 = lout[4]*1.0/u3p
				theta3 = u3p*1.0/(u3p+o3p)
			else:
				psi3 = "NA"
				if o3p > 0:
					theta3 = 0
				else:
					theta3 = "NA"

			lout += [u5p, u3p, o5p, o3p, psi5, psi3, theta5, theta3]
			print "\t".join(map(str, lout))

	else:

		# time to print out the junctions
		for jid in sorted(djuncs.keys()):
			print "\t".join(map(str, djuncs[jid][0:5]))


	return 0



#==============================================================================
# defs
#==============================================================================


#--
# parse_sam_header_row
# parses a sam header line into a list
def parse_sam_header_row(sz):
	aln = sz.strip().split("\t")
	dh = {}
	x = aln[0]
	dh["type"] = x[1:len(x)]
	ll = []
	for i in range(1,len(aln)):
		tmp = aln[i].split(":")
		ll.append(tmp[1])

	dh["data"] = ll
	return(dh)


#--
# parse_bam_alignment
# parse a bam alignment into a dict
def parse_bam_alignment(sz):
	# vars
	aln = []
	fields = ["qname", "flag", "rname", "pos", "mapq", "cigar", 
		"rnext", "pnext", "tlen", "seq", "qual", "ext"]
	ext_idx = 11
	daln = {}

	# explode
	aln = sz.strip().split("\t")
	nfields = len(aln)

	for i in range(0, ext_idx):
		daln[fields[i]] = aln[i]

	daln["ext"] = list(aln[ext_idx:nfields])

	return(daln)


def aln_is_spliced(daln):
	rres = False
	rr = re.search("[0-9]+N", daln["cigar"])
	if rr:
		rres = True
	return rres

#--
# explode_cigar
# explode cigar string into operations and lengths
# note: uses re
def explode_cigar(cigar):
	
	# using findall we can blow this thing up all in one shot
	res = []
	if cigar != "*":
		res = re.findall("([0-9]+)([MIDNSHPX=]+)", cigar)
	
	return(res)

#--
# aln_get_left_right
# get the left and right most coordinates (adding soft clips back in)
# for the alignment. expect a dict input.
def aln_get_left_right(daln):

	res = [-1, -1]
	left0 = 0
	left = 0
	right = 0

	if daln['cigar'] != "*" and aln_is_aln(daln):
		left0 = int(daln['pos'])
		left = left0
		right = left0

		# adjust left for soft-clipping?
		rres = re.search("^([0-9]+)S", daln["cigar"])
		if rres:
			left -= int(rres.group(1))

		# get ops that change position
		rres = re.findall("([0-9]+)([MDN])", daln["cigar"])
		# find right
		for i in range(len(rres)):
			right += int(rres[i][0])

		# adjust right for soft-clipping?
		rres = re.search("([0-9]+)S$", daln["cigar"])
		if rres:
			right += int(rres.group(1))

		right -= 1

		res = [left, right]

	return res


#--
# aln_get_mapped_regions
# get the mapped boundaries of the alignment
def aln_get_mapped_regions(daln):

	mrout = []
	left0 = 0
	left = 0
	right = 0

	if daln['cigar'] != "*" and aln_is_aln(daln):
		left = int(daln['pos'])

		# get ops that change position
		rres = re.findall("([0-9]+)([MDN=])", daln["cigar"])
		# find "M" ops or "=" ops
		for i in range(len(rres)):
			if rres[i][1] == "M" or rres[i][1] == "=":
				# make a region
				right = left+int(rres[i][0])-1
				mrout.append([daln["rname"], left, right])

			# move left position
			left += int(rres[i][0])

	return mrout

#--
# get a dict of True/False for all of the possible flags in an 
# alignment.
def aln_check_flags(daln):
	rres = {}
	flag = int(daln["flag"])

	rres["pe"] = (flag & 0x1)
	rres["proper"] = (flag & 0x2)
	rres["unmapped"] = (flag & 0x4)
	rres["mate_unmapped"] = (flag & 0x8)
	rres["reversed"] = (flag & 0x10)
	rres["mate_reversed"] = (flag & 0x20)
	rres["first_mate"] = (flag & 0x40)
	rres["second_mate"] = (flag & 0x80)
	rres["secondary"] = (flag & 0x100)
	rres["filtered"] = (flag & 0x200)
	rres["duplicate"] = (flag & 0x400)
	rres["supaln"] = (flag & 0x800)

	return(rres)


#--
# aln_is_aln
# return true if the alignment is aligned
def aln_is_aln(daln):
	res = True
	if int(daln["flag"]) & 0x4:
		res = False
	return(res)


#--
# hash_region
# r is a list with [ref, start, end] for the region
# and bin is a binning integer for hashing.
# the region may fall into multiple bins. 
# returns a list of hashes
def hash_region(r, bin):

	bin = int(bin)
	hstart = int(r[1])/bin
	hend = int(r[2])/bin

	hh = []
	for i in range(hstart, hend+1):
		hh.append("{}:{}".format(r[0], i))

	return(hh)

#--
# region_overlap
# compares regions a and b to see if they overlap. returns 
# a list: [0/1, overlap length]
def region_overlap(a, b):
	astart = float(a[1])

	if len(a) < 3:
		aend = astart
	else:
		aend = float(a[2])

	bstart = float(b[1])

	if len(b) < 3:
		bend = bstart
	else:
		bend = float(b[2])

	ovl = 0
	olen = 0

	# check overlap
	if aend >= bstart and astart <= bend:
		ovl = 1
		olen1 = aend-bstart
		olen2 = bend-astart
		len1 = aend-astart
		len2 = bend-bstart
		# overlap length is the minimum of all of these lengths
		olen = min([olen1, olen2, len1, len2])

	return([ovl, olen])

def prime_overlap(a, b):
	# special version of region overlap to implement finding hits to overlaps of 5' 
	# and 3' regions
	astart = float(a[1])
	aend = float(a[2])

	# b is a single point
	bstart = float(b[1])

	ovl = [0, 0]

	if astart < (bstart-4) and aend > (bstart+4):
		ovl[0] = 1
		ovl[1] = min([bstart-astart, aend-bstart])

	return(ovl)


#--
# build_point_lookup_hash
# lpoints is a list of, at minimum, two element lists that shall be hashed
def build_point_lookup_hash(lpoints, bin):

	pid = ""
	bin = int(bin)
	didx = {}
	
	for i in range(len(lpoints)):
		hpos = int(lpoints[i][1])/bin

		pid = "{}:{}".format(lpoints[i][0], hpos)

		if pid not in didx:
			didx[pid] = []

		# insert into the hash
		didx[pid].append(list(lpoints[i]))


	return(didx)

#--
# find_point_hash_hits
# look up hits to points for a range, a
def find_point_hash_hits(a, h, bin):

	# hash a, a region
	hh = hash_region(a, bin)
	# check for hits
	hids = {}

	for i in range(len(hh)):
		if hh[i] in h:
			# maybe, get list of elements at this node
			lcand = h[hh[i]]
			for j in range(len(lcand)):
				# overlap?
				rres = prime_overlap(a, lcand[j])
				if rres[0]==1:
					# keep this hit
					hid = ":".join(lcand[j])
					hids[hid] = 0

	return(hids.keys())



def range_to_id(r):
	rid = "{}:{}-{}".format(r[0], r[1], r[2])
	return(rid)

def point_to_id(r):
	pid = "{}:{}".format(r[0], r[1])
	return(pid)


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
parser.add_argument('bam', type=str, help="BAM alignments")
parser.add_argument('-l', type=int, default=180, 
	help="Estimate fragment length (for single end. paired-end is inferred from alignments) [180]")
parser.add_argument('-p', '--psi-theta', action="store_const", const=True, default=False,
	help="Compute psi/theta values [off]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

