#!/usr/bin/env python
#
# gtf-alternative-exons.py
#
# Shawn Driscoll
# 20130124
#
# part of the alternative exon process.  this script is run first on a GTF
# followed by running make-alt-exon-bed.py.
# 

import sys, argparse, re
import subprocess as sp
from hashlib import md5

#
# main function
#
def main(args):

	#
	# variables
	#
	final_list = []
	ae_id = ""
	out_file = args.gtf + ".ae"
	
	# check for the input file
	try:
		fin = open(args.gtf, 'r')
	except IOError,e:
		sys.stderr.write("Unable to open input file ({:s})\n".format(args.gtf))
		return(1)

	# check for the locus tag in the GTF file
	ll = fin.readline().strip().split("\t")
	attr = parse_gtf_attr(ll[8])
	
	if "locus" not in attr:
		#
		# check for cufflinks
		#
		# sys.stderr.write("checking for bowtie...\n")
		fnull = open('/dev/null', 'w')
		res = sp.call("which gffread", shell=True, stderr=fnull, stdout=fnull)
		fnull.close()
		if res == 1:
			sys.stderr.write("Error: cufflinks/gffread is not in the system path. Can't add locus tag.\n")
			return(1)
		
		# no locus tag, add one in with gffread
		sys.stderr.write("Locus tag was not found in your GTF, adding tag with gffread...\n")
		cmd = "gffread --cluster-only -TF -o {:s}.locus {:s}".format(args.gtf, args.gtf)
		fnull = open("/dev/null", "w")
		p1 = sp.Popen(cmd.split(), stdout=fnull)
		p1.wait()
		fnull.close()
		gtf_file = args.gtf + ".locus"
	else:
		gtf_file = args.gtf

	gtf = parse_gtf(gtf_file, type=args.feature_type)
	
	final_list = []
	# iterate through sorted locus list
	for lid in sorted(gtf.keys()):
		if len(gtf[lid].keys()) > 1:
			# this locus has multiple isoforms
			final_list += find_alt_exons_in_locus(gtf[lid])
	
	header = [ "ae_id",
		"ae_genoName",
		"ae_genoStart",
		"ae_genoEnd",
		"ae_transcript_ids",
		"ae_gene_names",
		"ae_exon_index",
		"strand",
		"exon_type",
		"ae_left",
		"ae_right",
		"skip_genoName",
		"skip_genoStart",
		"skip_genoEnd",
		"skip_transcript_ids",
		"skip_gene_names",
		"skip_intron_index",
		"skip_left",
		"skip_right"]
	
	# open output file
	fout = open(out_file, "w")
	
	fout.write("\t".join(header) + "\n")
	
	for i in range(len(final_list)):
		# add unique AE identifier to each row
		ae_id = "AELOC_{:08d}".format(i)
		fout.write(ae_id + "\t" + "\t".join(map(str, final_list[i])) + "\n")
	
	fout.close()
	
	

def find_alt_exons_in_locus(rloc):
	# rloc is a dict of all of the transcripts in a locus
	
	# get unique exons and introns
	exons = unique_locus_exons(rloc)
	introns = unique_locus_introns(rloc)
	
	if len(exons) == 0 or len(introns) == 0:
		return([])
		
	ci_skip = []
	ci_inc = []
	ae_list = []
	
	# loop through exons and look for introns that completely skip them
	for i in range(len(exons)):
		ei = exons[i]
		ei[1] = int(ei[1])
		ei[2] = int(ei[2])
		
		j = 0
		while introns[j][1] < ei[1]:
			# since loop only iterates if the start of the intron is to the left of the
			# current exon we only need to check if the end of the intron is to the right
			# of the end of the exon
			if introns[j][2] > ei[2]:
				# we have a skip, folks!
				ci_skip.append(j)
				ci_inc.append(i)
			
			j += 1
			# if j is past the end of the intron list break out
			if j == len(introns):
				break

	if len(ci_skip) == 0 and len(ci_inc) == 0:
		# no AE at this locus
		return([])
	
	# now we have lists of the indices of the skipping introns and the skipped exons. we can 
	# look up their information from the gtf data as well as find the flanking exons for each
	
	# loop through ckip indices (should be same count as inc indices)	
	for i in range(len(ci_skip)):
		skip_left = ""
		skip_right = ""
		inc_left = ""
		inc_right = ""
		exon_type = "internal"
		feature_strand = ""
		final_out = []		
		
		# get an isoform from the skip event
		skip_tid = introns[ci_skip[i]][3].split(";")[0]
		skip_idx = int(introns[ci_skip[i]][5].split(";")[0])
		
		# grab the exons that flank the intron within the same transcript. it doesn't matter
		# how many transcripts the exons are in - it's better if they are from the same transcript as
		# the intron only
		temp = list(rloc[skip_tid]['features'][skip_idx])
		skip_attr = parse_gtf_attr(temp[8])
		skip_left = [temp[0], temp[3], temp[4], skip_attr['transcript_id'], skip_attr['gene_id'], temp[-2]]
		
		temp = list(rloc[skip_tid]['features'][skip_idx+1])
		skip_attr = parse_gtf_attr(temp[8])
		skip_right = [temp[0], temp[3], temp[4], skip_attr['transcript_id'], skip_attr['gene_id'], temp[-2]]
		
		
		# grab the strand of the transcripts while we are at it
		feature_strand = rloc[skip_tid]['features'][skip_idx][6]
		
		# decode the hashes by finding which unique exons have the same hash then copying
		# them to the skip_left and skip_right variables
#		for k in range(len(exons)):
#			if exons[k][-1] == skip_left:
#				skip_left = list(exons[k][0:6])
#			elif exons[k][-1] == skip_right:
#				skip_right = list(exons[k][0:6])
		
		# get an isoform that includes the AE
		inc_tid = exons[ci_inc[i]][3].split(";")[0]
		inc_idx = int(exons[ci_inc[i]][5].split(";")[0])
		
		# grap the exon hashes for the two flanking exons to the AE. since it is possible
		# that the exon is a 3' or 5' there may not always be a left and right exon.
		if inc_idx > 0:
			inc_left = rloc[inc_tid]['features'][inc_idx-1][-1]
		if inc_idx+1 < len(rloc[inc_tid]['features']):
			inc_right = rloc[inc_tid]['features'][inc_idx+1][-1]
			
		# decode the hashes
		for k in range(len(exons)):
			if exons[k][-1] == inc_left:
				inc_left = list(exons[k][0:6])
			elif exons[k][-1] == inc_right:
				inc_right = list(exons[k][0:6])
		
		# tag exon type
		if len(inc_left) == 0:
			# no left exon
			if feature_strand == "+":
				exon_type = "5p"
			else:
				exon_type = "3p"
		elif len(inc_right) == 0:
			# no right exon
			if feature_strand == "+":
				exon_type = "3p"
			else:
				exon_type = "5p"
		
		# start building output for this event
		
		# append AE and related infos
		final_out = list(exons[ci_inc[i]][0:6])
		final_out.append(feature_strand)
		final_out.append(exon_type)
		final_out.append("|".join(map(str, inc_left)))
		final_out.append("|".join(map(str, inc_right)))
		# append skip and related infos
		final_out += list(introns[ci_skip[i]][0:6])
		final_out.append("|".join(map(str, skip_left)))
		final_out.append("|".join(map(str, skip_right)))
	
		ae_list.append(final_out)
	
	return(ae_list)
	
	
def unique_locus_exons(dl):

	# dl is a dict of transcript ids

	eset = set([]) # empty set for tracking unique exons
	lelist = []
	efinal = {}
	
	# pool exons
	for tid in dl.keys():
		if len(dl[tid]['features']) > 1:
			# only add exons from multiple exon isoforms
			temp = []
			temp = list(dl[tid]['features'])
			lelist += temp

	if len(lelist) < 3:
		# not enough exons to have differential splicing!
		return([])

	for i in range(len(lelist)):
		# update set with exon hash
		eset.update([lelist[i][-1]])

	# list of unique hashes
	ulist = list(eset)

	#
	# build a new hash of unique exons with some collapsed annotations
	#
	
	for i in range(len(lelist)):
		ehash = lelist[i][-1]
		if ehash in ulist:
			# this is an exon of interest
			attr = parse_gtf_attr(lelist[i][8])
			if ehash not in efinal:
				# grab the exon's geno info as well as the transcript id, gene id and exon index within the transcript
				efinal[ehash] = [lelist[i][0], lelist[i][3], lelist[i][4], attr['transcript_id'], attr['gene_id'], str(lelist[i][-2]), ehash]
			else:
				# append transcript id and gene id and exon index to the existing data for this exon
				efinal[ehash][3] += ";" + attr['transcript_id']
				efinal[ehash][4] += ";" + attr['gene_id']
				efinal[ehash][5] += ";" + str(lelist[i][-2])

	# transfer to a list and sort them by start position
	efinal_list = []
	
	for ehash in efinal.keys():
		efinal_list.append(list(efinal[ehash]))
	
	# sort list by start position of each exon row
	efinal_list.sort(key=lambda x: int(x[1]))

	# done!
	return(efinal_list)

#
# unique_locus_introns
# This function takes a locus and returns a list of all unique introns within the locus.
# single exon transcripts are ignored.
def unique_locus_introns(dl):

	# dl is a dict of transcript ids

	eset = set([])
	lelist = []
	efinal = {}
	efinal_list = []


	# generate and pool introns
	
	# iterate through transcripts
	for tid in dl.keys():
		# get exon list for transcript
		elist = dl[tid]['features']
		n = len(elist)
		if n > 1:
			# transcript has at least two exons so we can extract at least one intron
			for i in range(n-1):
				temp = list(elist[i])
				# end of this exon is start of intron
				temp[3] = int(temp[4])
				# start of next exon is end of intron
				temp[4] = int(elist[i+1][3])
				# set intron index
				temp[-2] = i
				lelist.append(temp)
			
	if len(lelist) < 1:
		# no introns = who cares, return empty list
		return([])

	for i in range(len(lelist)):
		#
		# make a hash for each intron from the start and end coordinates of the intron
		ehash = md5(str(lelist[i][3]) + str(lelist[i][4])).hexdigest()
		# set hash
		lelist[i][-1] = ehash

	#
	# build a new hash of introns with some collapsed annotations
	#

	for i in range(len(lelist)):
		ehash = lelist[i][-1]
		attr = parse_gtf_attr(lelist[i][8])
		if ehash not in efinal:
			# grab the introns's geno info as well as the transcript id, gene id and intron index within the transcript
			efinal[ehash] = [lelist[i][0], lelist[i][3], lelist[i][4], attr['transcript_id'], attr['gene_id'], str(lelist[i][-2])]
		else:
			# append transcript id and gene id and exon index
			efinal[ehash][3] += ";" + attr['transcript_id']
			efinal[ehash][4] += ";" + attr['gene_id']
			efinal[ehash][5] += ";" + str(lelist[i][-2])

	#
	# transfer to a list and sort them by start position
	#
	for ehash in efinal.keys():
		efinal_list.append(list(efinal[ehash]))
	
	# sort
	efinal_list.sort(key=lambda x: int(x[1]))

	# return list
	return(efinal_list)

#
# parse_gtf
# This function parses the GTF by locus and by transcript within each locus
def parse_gtf(gtf, type="exon"):

	# open the file
	fin = open(gtf, "r")

	dgtf = {}

	for szl in fin:

		ll = szl.strip().split("\t")
		
		if ll[2] == type:
			attr = parse_gtf_attr(ll[8])
	
			lid = attr['locus']
			tid = attr['transcript_id']
	
			# if locus id isn't present, add it
			if lid not in dgtf:
				dgtf[lid] = {}
	
	
			# if transcript id isn't present in the locus, add it
			if tid not in dgtf[lid]:
				dgtf[lid][tid] = {}
				dgtf[lid][tid]['features'] = []
	
			# append feature row list
			dgtf[lid][tid]['features'].append(ll)
	
	# close the GTF file
	fin.close()
	
	# make sure the exons in each transcript set are sorted by start position
	for lid in dgtf.keys():
		for tid in dgtf[lid].keys():
			# tid is current transcript id
			
			# sort the exon list within the transcript
			dgtf[lid][tid]['features'].sort(key=lambda x: int(x[3]))
		
			# now we can append the correct exon indices to the exon rows and we'll append a hash
			# value that'll be used later on to track exons
			for i in range(len(dgtf[lid][tid]['features'])):
				dgtf[lid][tid]['features'][i].append(i)
				ehash = md5(str(dgtf[lid][tid]['features'][i][3]) + str(dgtf[lid][tid]['features'][i][4])).hexdigest()
				dgtf[lid][tid]['features'][i].append(ehash)
		

	return(dgtf)


# 
# parse_gtf_attr
def parse_gtf_attr(field):

	attrs = {}

	# split into fields
	l1 = field.split(";")

	# split each field into key and value
	l2 = []
	for i in range(len(l1)):
		l2.append(l1[i].split("\""))


	for i in range(len(l2)):
		if len(l2[i]) > 1:
			attrs[l2[i][0].strip()] = l2[i][1].strip()

	return(attrs)


#
# main entry point
#

parser = argparse.ArgumentParser(description="Generates an alternative exon table from a GTF file. The output will be name of your GTF file with .ae appended to it.")
parser.add_argument('gtf', type=str, help="GTF file")
# options
parser.add_argument('-t', type=str, dest='feature_type', default="exon", help="Feature type (default: exon)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
	


