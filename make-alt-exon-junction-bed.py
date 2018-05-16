#!/usr/bin/env python
#
# make-alt-exon-bed.py
#
# Shawn Driscoll
# 20130124
#
# part of the alternative exon process.  this script looks for the .gtf.ae
# file generated with gtf-alternative-exons.py. this file you want to 
# quantify against transcriptome aligned reads.
# 

import sys, argparse, re
import subprocess as sp
from hashlib import md5

#
# main function
#
def main(args):

	# variables
	# fnull = open("/dev/null", "w")

	out_file = args.gtf + ".ae.junc.bed"
	alt_out_file = args.gtf + ".ae.exons.bed"

	# check the gtf file
	try:
		fin = open(args.gtf)
	except IOError,e:
		sys.stderr.write("Error: gtf file doesn't exist!\n")
		return(1)
	
	fin.close()
	
	# check for the AE file
	try:
		fin = open(args.gtf + ".ae")
	except IOError,e:
		sys.stderr.write("Can't find AE file that goes with your GTF file, I'll generate one now...\n")
		cmd = "gtf-alternative-exons {:s}".format(args.gtf)
		p1 = sp.Popen(cmd.split())
		p1.wait()
	else:
		fin.close()

	if args.exon_mode:
		# we only need to parse the AE file and generate a bed format for the AE's
		fin = open(args.gtf + ".ae", "r")
		fout = open(alt_out_file, "w")
		# skip header
		szl = fin.readline()
		for szl in fin:
#			print szl.strip()
			ll = szl.strip().split("\t")
			l_out = [ 
				ll[1], str(int(ll[2])-1), ll[3],
				"|".join(map(str,[ ll[0], ll[4], ll[5], ll[6], ll[7], ll[8] ])) ]
			fout.write("\t".join(l_out) + "\n")
			
		fout.close()
	
	else:

		# parse entire GTF into dict structure by transcript. this reference will be used
		# to help translate the genome cooridnates in the AE file into transcriptome coordinates
		# for hit counting from transcriptome alignments
#		gtf = parse_gtf(args.gtf, type=args.feature_type)
		
		# open and parse the AE file
		fin = open(args.gtf + ".ae", "r")
		fout = open(out_file, "w")
		# skip header
		szl = fin.readline()
		
		for szl in fin:
			ll = szl.strip().split("\t")
			
			# deal with the alternative exon
			
			# translate the coordinates for the exon and its flanking exons
			exon = ll[1:7]			
			exon_tid = exon[3].split(";")[0]
			exon_idx = int(exon[5].split(";")[0])
			# fetch this exon from the gtf
#			exon_t = gtf[exon_tid]['trans'][exon_idx]

			exon_left = len(ll[9].strip()) > 0
			exon_right = len(ll[10].strip()) > 0
			
#			info = [ exon_tid, str(exon_idx), ll[8] ]
			
			# write skip line
			skip = [ll[11], int(ll[12])+1, int(ll[13])-1]
			skip.append("%s.skip|%s|%s|intron" % (ll[0], ll[14], ll[16]))
			fout.write("\t".join(map(str, skip)) + "\n")
			
			inc_left = []
			inc_right = []
			if exon_left:
				temp = ll[9].split("|")				
				inc_left = [temp[0], int(temp[2])+1, int(ll[2])-1]
				inc_left.append("%s.inc1|%s|%s|%s" % (ll[0], ll[4], ll[6], ll[8]))
				fout.write("\t".join(map(str, inc_left)) + "\n")
				
			if exon_right:
				temp = ll[10].split("|")				
				inc_left = [temp[0], int(ll[3])+1, int(temp[1])-1]
				inc_left.append("%s.inc2|%s|%s|%s" % (ll[0], ll[4], ll[6], ll[8]))
				fout.write("\t".join(map(str, inc_left)) + "\n")
			
		
		fin.close()
		fout.close()
		
		# we need the bed file to be sorted by columns 1 and 2 so we'll do that 
		# with the system. we can also collapse the list based on matching 
		# positions
#		fout = open(out_file + ".sort", "w")		
#		p1 = sp.Popen("sort -k1,1 -k2,2n {:s}".format(out_file).split(), stdout=sp.PIPE)
#		p2 = sp.Popen("groupBy -i stdin -g 1,2,3 -c 4 -o collapse".split(), stdin=p1.stdout, stdout=fout)
#		p2.wait()
#		fout.close()
		# move the sorted file over the original
#		p1 = sp.Popen("mv {:s} {:s}".format(out_file + ".sort", out_file).split())
#		p1.wait()
		# done!

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
	
			tid = attr['transcript_id']
	
	
			# if transcript id isn't present in the locus, add it
			if tid not in dgtf:
				dgtf[tid] = {}
				dgtf[tid]['features'] = []
	
			# append feature row list to transcript within locus
			dgtf[tid]['features'].append(ll)
			
	fin.close()

	# make sure that for each transcript list the exon features are sorted by start position. 
	# then we'll add the exon index to each one
	
	for tid in dgtf.keys():
		# tid is current transcript id
		
		# sort the exon list within the transcript
		dgtf[tid]['features'].sort(key=lambda x: int(x[3]))
		
		# now for each transcript we want to generate transcriptome coordinates
		elist = dgtf[tid]['features']
		dgtf[tid]['trans'] = []
		
		elens = []
		
		for i in range(len(elist)):
			elens.append(int(elist[i][4]) - int(elist[i][3]))
			dgtf[tid]['trans'].append(list(elist[i]))
		
		last = 0
		for i in range(len(elist)):
			dgtf[tid]['trans'][i][3] = last+1
			dgtf[tid]['trans'][i][4] = dgtf[tid]['trans'][i][3]+elens[i]
			last = dgtf[tid]['trans'][i][4]
		

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

parser = argparse.ArgumentParser(description="Uses an alt-exon table (from gtf-alternative-exons.py) to generate a bed format file that can be used to quantify skip/inclusion hits from transcriptome alignments")
parser.add_argument('gtf', type=str, help="GTF file for quantification")

# options
parser.add_argument('-t', type=str, dest='feature_type', default="exon", help="Feature type (default: exon)")
parser.add_argument("--exon-mode", dest='exon_mode', action="store_const", const=True, default=False, help="Exons mode produces a bed file of the alternative exons only. (default: off)")
parser.add_argument("-w", dest="overlap_window", type=int, default=5, help="Size of overlap window for overlapping read detection at junctions. Indicates the number of bases past the junction in each direction. (default: 5)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
	


