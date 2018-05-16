#!/usr/bin/env python
#
# bwa-express.py
#
# Shawn Driscoll
# 20121109
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# align for eXpress expression analysis
#

import sys,argparse
import subprocess as sp

#
# main
#

# example: bwa aln -n 0.5 -l 25 -k 2 -t 8 -R 200 /path/to/index.fa left.reads.fastq > left.sai
#          bwa aln -n 0.5 -l 25 -k 2 -t 8 -R 200 /path/to/index.fa right.reads.fastq > right.sai
#          bwa sampe -a 1000 -n 200 -N 0 /ath/to/index.fa left.sai right.sai left.reads.fastq right.reads.fastq > aligned.sam
# for single end:
#          bwa samse -n 200 /path/to/index.fa left.sai left.reads.fastq > aligned.sam

def main(args):

	cmd_express = None

	#
	# check for commands
	#
	sys.stderr.write("checking for bwa...\n")
	res = sp.call("which bwa",shell=True,stderr=sp.STDOUT)
	if res == 1:
		sys.stderr.write("Error: bowtie is not in the system path\n")
		return 1

	if not args.pipe_to_express:
		sys.stderr.write("checking for samtools...\n")
		res = sp.call("which samtools",shell=True,stderr=sp.STDOUT)
		if res == 1:
			sys.stderr.write("Error: samtools is not in the system path\n")
			return 1

	sys.stderr.write("checking for eXpress...\n")
	res = sp.call("which express",shell=True,stderr=sp.STDOUT)
	if res == 1:
		sys.stderr.write("Error: express is not in the system path\n")
		return 1

	#
	# check for reads
	#
	if args.paired and len(args.reads) < 2:
		sys.stderr.write("Error: not enough read files specified for paired alignments!\n")
		return 1

	#
	# figure out if the reads are gzip compressed or not
	#
	res = sp.check_output("file {:s}".format(args.reads[0]),shell=True)
	if res.find("gzip") >= 0:
		#
		# reads are gzipped so they'll need to be gunzipped when bowtie is run
		#
		for i in range(len(args.reads)):
			args.reads[i] = "<(gunzip -c " + args.reads[i] + ")"


	#
	# start building command strings
	#

	if args.paired:
		aln_stub = "bwa aln -n 0.5 -l 25 -k 2 -t {:d} -R 200 {:s} {:s} > {:s}"
		aln_cmd = aln_stub.format(args.num_threads,args.index,args.reads[0],"left_tmp.sai")
		aln_cmd = aln_cmd + "\n" + aln_stub.format(args.num_threads,args.index,args.reads[1],"right_tmp.sai")
		aln_cmd = aln_cmd + "\n" + "bwa sampe -a 1000 -n 200 -N 0 {:s} {:s} {:s} {:s} {:s}".format(args.index,"left_tmp.sai","right_tmp.sai",args.reads[0],args.reads[1])
	else:
		aln_stub = "bwa aln -n 0.5 -l 25 -k 2 -t {:d} -R 200 {:s} {:s} > {:s}"
		aln_cmd = aln_stub.format(args.num_threads,args.index,args.reads[0],"aln_tmp.sai")
		aln_cmd = aln_cmd + "\n" + "bwa samse -n 200 {:s} {:s} {:s}".format(args.index,"aln_tmp.sai",args.reads[0])

	if args.pipe_to_express:
		#
		# pipe alignments straight to eXpress
		#
		aln_cmd = aln_cmd + " | express -m {:d} -s {:d} {:s}".format(args.fragment_mean,args.fragment_sd,args.index)
	else:
		aln_cmd = aln_cmd + " | samtools view -bS -F 0x04 -o {:s}.bam -".format(args.stub)

		# skip expression?
		if not args.no_expr:
			aln_cmd = aln_cmd + "\nexpress -m {:d} -s {:d} {:s} {:s}.bam".format(args.fragment_mean,args.fragment_sd,args.index,args.stub)

		# kill bam file afterwards?
		if args.no_bam:
			aln_cmd = aln_cmd + "\nrm {:s}.bam".format(args.stub)

	#
	# remove .sai files
	#
	if args.paired:
		aln_cmd = aln_cmd + "\nrm left_tmp.sai right_tmp.sai"
	else:
		aln_cmd = aln_cmd + "\nrm aln_tmp.sai"

	if args.echo_only:
		print aln_cmd
	else:
		#
		# write commands out to a temporary script file
		#
		fp = open("py_temp.sh","w")
		fp.write("#!/bin/bash\n")
		fp.write(aln_cmd + "\n")
		fp.close()

		#
		# chmod and run the temporary script
		#
		sp.call("chmod 755 py_temp.sh && ./py_temp.sh",shell=True)

		#
		# see if stuff finished
		#
		try:
			fp = open("results.xprs",'r')
		except IOError,e:
			sys.stderr.write("it doesn't look like eXpress did its thing\n")
			return 1

		fp.close()
		sp.call("mv results.xprs {:s}.results.xprs".format(args.stub),shell=True)
		sp.call("mv params.xprs {:s}.params.xprs".format(args.stub),shell=True)
		sp.call("rm py_temp.sh",shell=True)

	return 0

#==============================================================================
# entry point
#==============================================================================

# modes list

# bowtie-express
# bowtie2-express
# bwa-express
# tophat-cufflinks

parser = argparse.ArgumentParser(description="g3nom1c!")
parser.add_argument("reads",metavar="read_files",type=str,nargs="+",help="reads. if multiple read files they must be separated by a single space (send xargs to this command)")
parser.add_argument("-o",dest="stub",type=str,default="aligned",help="stub for output files (default: aligned)")
parser.add_argument("-x",dest="index",type=str,default="mm9_kg",help="bowtie index (default: mm9_kg)")
parser.add_argument("--paired",dest="paired",action="store_const",const=True,default=False,help="data is paired end")
parser.add_argument("-p",dest="num_threads",type=int,default=8,help="number of threads for alignment (default: 8)")
parser.add_argument("--no-expr",dest="no_expr",action="store_const",const=True,default=False,help="skip expression quantification")
parser.add_argument("--pipe",dest="pipe_to_express",action="store_const",const=True,default=False,help="pipe alignments directly to eXpress (no bam file retained)")
parser.add_argument("--no-bam",dest="no_bam",action="store_const",const=True,default=False,help="do not keep the bam alignments file")
parser.add_argument("-m",dest="fragment_mean",type=int,default=350,help="mean fragment length for eXpress (default: 350)")
parser.add_argument("-s",dest="fragment_sd",type=int,default=80,help="fragment length standard deviation (default: 80)")
parser.add_argument("--echo-only",dest="echo_only",action="store_const",const=True,default=False,help="echo command string to stdout, don't execute it. (default: False)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
