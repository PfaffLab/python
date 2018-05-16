#!/usr/bin/env python
#
# transcriptome-counts.py
#
# Shawn Driscoll
# 20121215
#
# Aligns reads to a transcriptome (pre-built index required) and sums counts 
# into loci based on a GTF that's been parsed into loci (by gffread)
#
# 20130329
# modified the file to run bowtie2 instead of bowtie
#

import sys, argparse, re
import subprocess as sp

#
# main function
#
def main(args):

	# variables
	num_files = len(args.read_files)
#	bowtie_options = "-p {:d} -n {:d} -e {:d} -M {:d} --best -t -S ".format(args.num_threads, args.seed_mismatches, args.max_qual_mismatch, args.max_multi_report)
	bowtie_options = "-p {:d} -N {:d} ".format(args.num_threads, args.seed_mismatches)

	ltable = {}
	rtable = {}
	has_locus = False
	has_gene_name = False
	
	#
	# check all files
	#
	try:
		sys.stderr.write("> checking gtf...")
		fin = open(args.gtf, "r")
	except IOError,e:
		sys.stderr.write(e.strerror + "\n")
		return e.errno
	else:
		sys.stderr.write("found!\n")
		fin.close()
	
	try:
		sys.stderr.write("> checking bowtie reference...")
		#fin = open(args.tref + ".1.ebwt", "r")
		fin = open(args.tref + ".1.bt2", "r")
	except IOError,e:
		sys.stderr.write(e.strerror + "\n")
		return e.errno
	else:
		sys.stderr.write("found!\n")
		fin.close()
	
	for i in range(num_files):
		try:
			sys.stderr.write("> checking read file {:s}...".format(args.read_files[i]))
			fin = open(args.read_files[i],"r")
		except IOError,e:
			sys.stderr.write(e.strerror + "\n")
			return e.errno
		else:
			sys.stderr.write("found!\n")
			fin.close()

	#
	# check for locus map file for GTF
	#
#	try:
#		fin = open(args.gtf + ".lt", "r")
#	except IOError:
#		sys.stderr.write("> cannot find locus info table for your GTF, building one now...\n")
#		cmd = "gtf-loci --no-group {:s} > {:s}.lt".format(args.gtf, args.gtf)
#		sp.call("bash -c '{:s}'".format(cmd), shell=True)
#		fin = open(args.gtf + ".lt", "r")

	fin = open(args.gtf, "r")
	for szl in fin:
		ll = szl.strip().split("\t")
		if ll[2] == "exon":
			m = re.search('locus "[^"]+"', szl)
			if m:
				has_locus = True
			else:
				has_locus = False

			break

	fin.close()

	if not has_locus:
		sys.stderr.write("> your gtf doesn't have the locus tag. adding it now...\n")
		cmd = "locusify-gtf {:s}".format(args.gtf)
		sp.call("bash -c '{:s}'".format(cmd), shell=True)		
	
	sys.stderr.write("> hashing locus table...\n")
	fin = open(args.gtf, "r")
	for szl in fin:
		szl = szl.strip()
		if len(szl) > 0:
			ll = szl.split("\t")
			if ll[2] == "exon":
				attr = parse_gtf_attr(ll[8])

				if attr['transcript_id'] not in ltable:
					ltable[attr['transcript_id']] = [ attr['locus'], attr['gene_id'], ""]
					if "gene_name" in attr:
						has_gene_name = True
						ltable[attr['transcript_id']][2] = attr['gene_name']

	fin.close()

	#
	# loop through files and do the job
	#
	for i in range(num_files):
		#
		# build stub for output
		#
		ltmp = args.read_files[i].split("/")
		ltmp = ltmp[len(ltmp)-1].split(".")
		stub = ltmp[0]

		# check for gzipped reads		
#		res = sp.check_output("file {:s}".format(args.read_files[i]),shell=True)
#		m = re.search("gzip",res)
		
#		if m:
#			rfile = "<(gunzip -c {:s})".format(args.read_files[i])
#		else:
#			rfile = args.read_files[i]
		
		rfile = args.read_files[i]
		
		# run bowtie
		sys.stderr.write("> aligning {:s} with bowtie2...\n".format(args.read_files[i]))
#		cmd = "bowtie {:s} {:s} {:s} 2> {:s}.bwtlog | samtools view -bS -F 0x4 -o temp.bam -".format(bowtie_options, args.tref, rfile, stub)
		cmd = "bowtie2 {:s} -x {:s} -U {:s} 2> /dev/null | samtools view -bS -F 0x4 -o temp.bam - 2> /dev/null".format(bowtie_options, args.tref, rfile)
		sp.call("bash -c '{:s}'".format(cmd), shell=True)
		
		sys.stderr.write("> sorting alignments...\n")
		sp.call("samtools sort temp.bam {:s}".format(stub), shell=True)
		sp.call("samtools index {:s}.bam".format(stub), shell=True)
		# sp.call("bash -c 'samtools idxstats {:s}.bam > temp.counts'".format(stub), shell=True)
		
		
		sys.stderr.write("> counting hits and merging locus information...\n")
		rtable = {}
		p1 = sp.Popen("samtools idxstats {:s}.bam".format(stub).split(), stdout=sp.PIPE)
		
		# read in hits, fetch locus info to sum hits into loci
		
		for szl in p1.stdout:
			szl = szl.strip()
			if len(szl) > 0:
				arl = szl.split("\t")
				# look up feature in hash
				if arl[0] in ltable:
					tid = arl[0]
					lid = ltable[tid][0]
					# update final hash
					if lid not in rtable:
						rtable[lid] = {}
						rtable[lid]["tids"] = set([])
						rtable[lid]["gids"] = set([])
						if has_gene_name:
							rtable[lid]["gns"] = set([])
						rtable[lid]['counts'] = 0
					
					rtable[lid]['tids'].update([tid])
					rtable[lid]['gids'].update([ltable[tid][1]])
					if has_gene_name:
						rtable[lid]["gns"].update([ltable[tid][2]])
					rtable[lid]['counts'] += int(arl[2])
		
		#
		# write final locus hits out to file
		#
		
		sys.stderr.write("> writing final counts...\n")
		fout = open(stub + ".counts", "w")
		
		# write header
		if has_gene_name:
			fout.write("locus_id\tgene_names\tgene_ids\ttranscript_names\tcounts\n")
		else:
			fout.write("locus_id\tgene_names\ttranscript_names\tcounts\n")
		
		for lid in sorted(rtable.keys()):
			fout.write(lid + "\t")
			if has_gene_name:
				fout.write(",".join(sorted(list(rtable[lid]['gns']))) + "\t")
			fout.write(",".join(sorted(list(rtable[lid]['gids']))) + "\t")
			fout.write(",".join(sorted(list(rtable[lid]['tids']))) + "\t")
			fout.write(str(rtable[lid]['counts']) + "\n")
		
		fout.close()
		
		# remove temp files
		sp.call("rm temp.bam", shell=True)
		
		if not args.keep_bam:
			sp.call("rm {:s}.bam".format(stub), shell=True)
		
		sys.stderr.write("> done!\n")

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

parser = argparse.ArgumentParser(description="Aligns reads to a transcriptome. And quantifies hits per locus.")
parser.add_argument('tref', type=str, help='transcriptome reference (bowtie1)')
parser.add_argument('gtf', type=str, help="GTF to match the transcriptome being alinged to")
parser.add_argument('read_files',type=str,metavar='reads',nargs='+',help='FASTQ reads')
# options
parser.add_argument('-p', type=int, dest='num_threads', default=1, help="Number of threads to use during alignment (default: 1)")
parser.add_argument('-n', type=int, dest='seed_mismatches', default=1, help="Number of mismatches allowed in the seed (0-2) (default: 1)")
parser.add_argument('-e', type=int, dest='max_qual_mismatch', default=80, help="-e option in Bowtie (default: 80)")
parser.add_argument("-M", type=int, dest="max_multi_report", default=100, help="-M option in Bowtie (default: 100)")
parser.add_argument("--keep-bam", dest="keep_bam", action="store_const", const=True, default=False, help="Keep alignments (default: False)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
	
