#!/usr/bin/env python
#
# express-annotate-xprs.py
#
# Shawn Driscoll
# 20130128
#
# Pulls gene names and locus indexes from gtf information. Generates an isoform
# level file will all available columns and a gene locus level file with 
# only the columns that can be merged per locus (counts, fpkm).
#

import sys, argparse
import subprocess as sp

#
# MAIN
#
def main(args):

	# first we should check for the ".lt" output of gtf-loci for annotation
	try:
		fin = open(args.gtf + ".lt")
	except:
		sys.stderr.write("Unable to locate locus info table. I'll build one now.\n")
		fout = open(args.gtf + ".lt", "w")
		cmd = "gtf-loci --no-group {:s}".format(args.gtf)
		p1 = sp.Popen(cmd.split(), stdout=fout)
		p1.wait()
		fout.close()
	else:
		fin.close()

	# parse in the locus info table by transcript id
	transcript_info = {}
	locus_info = {}
	fin = open(args.gtf + ".lt", "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		transcript_info[ll[1]] = [ll[0], ll[2]]

		if ll[0] not in locus_info:
			# each locus will have a list of unique transcript ids and unique gene names
			locus_info[ll[0]] = [set([]), set([])]

		# update locus info
		locus_info[ll[0]][0].update([ll[1]])
		locus_info[ll[0]][1].update([ll[2]])

	fin.close()

	# loop through files and annotate/collapse, etc
	for i in range(len(args.xprs_files)):
		fin = open(args.xprs_files[i], "r")

		locus_hits = {}

		# read in header
		ll = fin.readline().strip().split("\t")
		# modify the header
		new_head = ["locus_id", "gene_name"] + list(ll[1:])
		locus_head = ["locus_id", "transcript_ids", "gene_name", "counts", "fpkm"]

		fstub = args.xprs_files[i].split(".")[0]
		# we can write the isoform level file out in this loop
		fout = open(fstub + ".isoforms.xprs", "w")
		# print the header
		fout.write("\t".join(new_head) + "\n")

		# loop through file and build outputs
		for szl in fin:
			ll = szl.strip().split("\t")

			# deal with locus level info
			locus_id = transcript_info[ll[1]][0]

			if locus_id not in locus_hits:
				# eff_counts and fpkm
				locus_hits[locus_id] = [0.0, 0.0]

#			locus_hits[locus_id][0] += float(ll[6])
			locus_hits[locus_id][0] += float(ll[7])
			locus_hits[locus_id][1] += float(ll[10])

			# write the isoform level row
			ll_new = [locus_id, transcript_info[ll[1]][1]] + list(ll[1:])
			fout.write("\t".join(ll_new) + "\n")

		fout.close()

		# now we can write the gene locus level summary
		fout = open(fstub + ".genes.xprs", "w")
		# write the header
		fout.write("\t".join(locus_head) + "\n")
		for lid in sorted(locus_hits.keys()):
			lout = [ lid, 
				",".join(list(locus_info[lid][0])), 
				",".join(list(locus_info[lid][1])), 
				str(locus_hits[lid][0]), str(locus_hits[lid][1]) ]

			fout.write("\t".join(lout) + "\n")

		fout.close()





#
# main entry point
#

parser = argparse.ArgumentParser(description="Pulls gene names and locus indexes from gtf information. Generates an isoform level file will all available columns and a gene locus level file with only the columns that can be merged per locus (counts, fpkm).")
parser.add_argument('gtf', type=str, help="GTF annotation corresponding to quantified data.")
parser.add_argument('xprs_files', type=str, metavar='xprs', nargs='+', help='xprs quantification file(s)')
# options

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
	