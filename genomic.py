#!/usr/bin/env python
#
# genomic.py
#
# Shawn Driscoll
# 20121109
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# one portal for all pipelines
#

import sys,argparse
import subprocess as sp

#
# main
#


def main(args):

	#cmd = ["bowtie","-p 8", "-n 2", "-a", "-m 200", "-t", "-S"]
	#flog = open("bowtie.log","w")

	#sp.call(cmd + [args.index,args.se_reads,"|","samtools","view","-bS","-o aligned.bam","-"],stderr=flog)

	#sp.call("bowtie -p 8 -n 2 -a -m 200 -t -S " + args.index + " " + args.se_reads + " | samtools view -bS -F 0x04 -o aligned.bam -",shell=True)

	#flog.close()

	if args.mode == "bowtie-express":


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
parser.add_argument('mode',type=str,help='operation mode (see list below)')
parser.add_argument('--source-dir',dest="source_dir",default=".",help="source folder for input files")
parser.add_argument('--reads',metavar="reads",type=str,nargs="+",help="reads in FASTQ for alignment modes")
parser.add_argument
parser.add_argument('-x',dest="index",type=str,help="bowtie index")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
