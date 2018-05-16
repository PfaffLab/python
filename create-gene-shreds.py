#!/usr/bin/python
#==============================================================================
# create-gene-shreds.py
#
# Shawn Driscoll
# 20180201
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
from Basics import messages as ms
from os import system, unlink
import hashlib
from multiprocessing import Process, JoinableQueue, current_process, Queue
import pysam as ps
from GenomeJunk import PysamTools as ps_tools
import numpy as np

# from igraph import *
# from subprocess import Popen
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

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	# check input files
	if not isfile(args.fasta):
		ms.error_message("Input FASTA file does not exist")
		return 1
	
	if not isfile(args.ref):
		ms.error_message("Input refFlat annotation does not exist")
		return 1
	
	if args.o is None:
		args.o = "{}.shreds.{}.fasta".format(drop_file_ext(args.fasta), args.l)
	else:
		tmp = "{}.shreds.{}.fasta".format(args.o, args.l)
		args.o = tmp
	
	if not args.just_alignment and not args.just_quantification:
	
		if isfile(args.o):
			unlink(args.o)
		
		if isfile("{}.gz".format(args.o)):
			unlink("{}.gz".format(args.o))
	
	if isfile("shred.log"):
		unlink("shred.log")
	
	if isfile("dedupe.log"):
		unlink("dedupe.log")
	
	if args.quantify:
		args.bbmap = True
	
	
	# good to go!
	
	if args.s:
		rres = core2(args)
	else:
		rres = core(args)

	return rres

def core2(args):
	##
	## in this mode we only have to align the input file to the specified index
	
	bam_out = re.sub("\.fasta$", ".bam", args.fasta)
	rres = bbmap(args.fasta, args.ref, bam_out)
	
	return 0
	

def core(args):
	
	annot = {}
	tid2gname = {}
	
	stub = hashlib.md5(args.fasta).hexdigest()
	
	fixed_case = "{}.ref.fa".format(stub)
	#gene_seq = "{}.gene.fa".format(stub)
	gene_shred = "{}.shred.fa".format(stub)
	gene_final = "{}.final.fa".format(stub)
	sub_chars = re.escape("[]{}\|/?!@#$%^&*()+=.") + "\s"
	
	fset = [fixed_case, "{}.fai".format(fixed_case), gene_shred, gene_final]

	##
	# parse refflat
	annot, tid2gname = load_refflat(args.ref)
	
	##
	# change case of input fasta
	ms.message("Converting all reference bases to uppercase")
	rres = fasta_fix_case(args.fasta, fixed_case)

	if args.just_alignment:
		args.bbmap = True
	
	if args.just_quantification:
		args.quantify = True
	
	if not args.just_alignment and not args.just_quantification:
		
		##
		# creat worker process to deal with all of the shredding
		tasks = JoinableQueue()
		results = Queue()
		p = None
		pool = []
		
		for i in range(args.p):
			p = Process(target=worker, args=(tasks, results, ))	
			p.daemon = True
			p.start()	
			pool.append(p)
			
		##
		# get to work!
		i = 0
		n = len(annot.keys())
		ms.message("Starting main loop")
		for gid in annot.keys():
			i += 1
			if (i % 3) == 0:
				ms.progress_message("processing {}. {} of {}".format(gid, i, n))
			
			gidHat = re.sub("[{}]".format(sub_chars), "", gid)
			
			# export all sequences belonging to this gene
			gene_seq = "{}.{}.fa".format(stub, gidHat)
			
			rres = samtools_faidx(fixed_case, gene_seq, annot[gid])
			tasks.put([stub, gidHat])
		
		for p in pool:
			tasks.put(None)
		
		ms.progress_message("Waiting for shredding to complete", last=True)
		
		tasks.join()
		for p in pool:
			p.join()
		
		ms.message("done")
		
		results.put(None)
		
		while True:
			fname = results.get()
			if fname is None:
				break
			
			if not isfile(fname):
				continue
			
			ms.message("Joining {}".format(fname))
			
			rres = cat_result(fname, args.o)
			unlink(fname)
		
	# done!
	
	bam_out = re.sub("\.fasta$", ".bam", args.o)
	
	if args.bbmap:
		if not isfile(args.o):
			if not isfile("{}.gz".format(args.o)):
				ms.error_message("Shredded reads file doest not exist, gzipped or not ({})".format(args.o))
				return 1
			else:
				# gunzip the reads file
				runcmd("gunzip {}.gz".format(args.o))
				
		ms.message("Aligning shreds back against the reference")
		
		rres = bbmap(args.o, fixed_case, bam_out, args.t)
	
	if args.quantify:
		if not isfile(bam_out):
			ms.error_message("Expected alignment file does not exist ({})".format(bam_out))
			return 1
			
		ms.message("Parsing alignments")
		pares = parse_alignments(bam_out, tid2gname)
		rres = process_pares(pares)
		
		tsvout = re.sub("\.bam$", ".tsv", bam_out)
		
		# output:
		# gene_name, total_reads, min_mapp, max_mapp, mean_mapp, most_similar, all_genes, all_gene_counts
		ms.message("Writing mappability report")
		with open(tsvout, "w") as fout:
			
			fout.write("#read_length={}\n".format(args.l))
			fout.write("\t".join(["#gene_name", "total_reads", "min_mapp", "max_mapp", "mean_mapp", "most_similar", "all_targets", "target_counts"]) + "\n")
			
			for gname in sorted(rres.keys()):
				
				total_reads = rres[gname]['total_reads']
				if total_reads==0:
					ms.warning_message("{} had zero reads".format(gname))
				min_mapp = 1
				max_mapp = 1
				mean_mapp = 1
				most_similar = "na"
				all_genes = "na"
				all_gene_counts = "na"
				
				if len(rres[gname]['target']) > 0:
					mapp = []
					for n in rres[gname]['target_count']:
						if total_reads > 0:
							mapp.append(1 - n*1.0/total_reads)
						else:
							mapp.append(0)
					
					min_mapp = min(mapp)
					max_mapp = max(mapp)
					mean_mapp = np.mean(mapp)
					
					for i in range(len(mapp)):
						if mapp[i]==min_mapp:
							most_similar = rres[gname]['target'][i]
					
					all_genes = ",".join(rres[gname]['target'])
					all_gene_counts = ",".join(map(str, rres[gname]['target_count']))
				
				lout = [gname, total_reads, min_mapp, max_mapp, mean_mapp, most_similar, all_genes, all_gene_counts]
				fout.write("\t".join(map(str, lout)) + "\n")
	
	if not args.just_quantification and not args.just_alignment:
		if args.z:
			if isfile("{}.gz".format(args.o)):
				unlink("{}.gz".format(args.o))
				
			cmd = "gzip {}".format(args.o)
			runcmd(cmd)
	
	#
	# clear out temp files
	for fname in fset:
		if isfile(fname):
			unlink(fname)

	
	return 0

def worker(qin, qout):
	name = current_process().name
	
	outfile = "{}.reads.fa".format(name)
	# create and or clear out the output file for this process.
	runcmd("cat /dev/null > {}".format(outfile))
	
	gene_shred = "{}.shred.fa".format(name)
	gene_final = "{}.final.fa".format(name)
	
	if isfile(outfile):
		unlink(outfile)
	
	while True:
		task = qin.get()
		if task is None:
			qin.task_done()
			break
		
		##
		# the task is a list with the stub and gene name
		stub = task[0]
		gname = task[1]
		
		gene_seq = "{}.{}.fa".format(stub, gname)

		if isfile(gene_seq):

			# shred the sequences
			rres = shred(gene_seq, gene_shred, args.l)
			# dedupe
			rres = dedupe(gene_shred, gene_final)
			# concat
			rres = cat_result(gene_final, outfile)
			
			unlink(gene_seq)
		else:
			ms.warning_message("{} is missing".format(gene_seq))
		
		qin.task_done()
	
	# send this thread's output file name back to main
	qout.put(outfile)
	
	for fname in [gene_shred, gene_final]:
		if isfile(fname):
			unlink(fname)
	
	return 0

def parse_alignments(fname, annot):
	
	# re.sub("\_[0-9]+\-[0-9]+$", "", sz)

	phash = {}

	qname = ""
	qname_last = ""
	rcount = 0

	with ps.AlignmentFile(fname) as fin:
		rnames = ps_tools.get_alignmentfile_rnames(fin)
		qbuff = set()
		
		for aln in fin:
			if aln.flag & 0x4:
				continue
			
			qname = aln.query_name
			if qname != qname_last:
				qbuff = list(qbuff)
				
				if len(qbuff) > 0:
					# deal with it
					
					rcount += 1
					if (rcount % 1000000) == 0:
						ms.progress_message("parsed {} reads".format(rcount))
					
					for k in qbuff:
						tmp = k.split("|")
						g_source = tmp[0]
						g_target = tmp[1]

						if g_source not in phash:
							phash[g_source] = {}
						
						if g_target not in phash[g_source]:
							phash[g_source][g_target] = 0
						
						phash[g_source][g_target] += 1
				
				qbuff = set()
							
			# get the target
			read_source = re.sub("\_[0-9]+\-[0-9]+$", "", aln.query_name)
			target = rnames[aln.reference_id]
			
			g_source = annot[read_source]
			g_target = annot[target]
			
			qbuff.add("{}|{}".format(g_source, g_target))
			qname_last = qname
			
	# out of loop. handle the final read
	
	qbuff = list(qbuff)	
	if len(qbuff) > 0:
		# deal with it

		for k in qbuff:
			tmp = k.split("|")
			g_source = tmp[0]
			g_target = tmp[1]

			if g_source not in phash:
				phash[g_source] = {}
			
			if g_target not in phash[g_source]:
				phash[g_source][g_target] = 0
			
			phash[g_source][g_target] += 1
	
	ms.progress_message("parsed {} reads [fin]".format(rcount), last=True)	
	
	return phash

def process_pares(d):
	##
	# build a report
	
	dreport = {}
	
	for gid in d.keys():
		if gid not in dreport:
			dreport[gid] = { "total_reads":0, "target":[], "target_count":[] }
			
		for tid in d[gid].keys():
			if gid == tid:
				# this is the alignments to itself. keep this count
				dreport[gid]["total_reads"] = d[gid][tid]
			else:
				dreport[gid]["target"].append(tid)
				dreport[gid]["target_count"].append(d[gid][tid])
				
	return dreport
		

def drop_file_ext(fname):
	tmp = fname.split(".")
	return ".".join(tmp[0:(len(tmp)-1)])

def cat_result(infile, target):
	cmd = "cat {} >> {}".format(infile, target)
	return runcmd(cmd)

##
# we only need a gene name to transcript name translation table
def load_refflat(fname):
	
	gene2tid = defaultdict(list)
	tid2gene = {}
	
	with open(fname, "r") as fin:
		
		for szl in fin:
			if szl[0] == "#":
				continue
			
			aln = szl.strip().split("\t")
			gene2tid[aln[0]].append(aln[1])
			tid2gene[aln[1]] = aln[0]
	
	return gene2tid, tid2gene

def fasta_fix_case(fname_in, fname_out):
	
	with open(fname_in, "r") as fin, open(fname_out, "w") as fout:
		for szl in fin:
			if szl[0]==">":
				fout.write(szl)
			else:
				fout.write(szl.strip().upper() + "\n")
	
	cmd = "samtools faidx {}".format(fname_out)
	runcmd(cmd)
	
	return 0

def samtools_faidx(ref, outfile, tlist):
	sz = " ".join(tlist)
	cmd = "samtools faidx {} {} > {}".format(ref, sz, outfile)
	runcmd(cmd)
	return 0
	
def dedupe(infile, outfile):
	cmd = "dedupe.sh in={} out={} overwrite threads=1 2>>dedupe.log".format(infile, outfile)
	runcmd(cmd)
	return 0

def bbmap(infile, ref, outfile, threads):
	cmd = "bbmap.sh in={} ref={} out={} ambig=all threads={} perfectmode trd tuc".format(infile, ref, outfile, threads)
	runcmd(cmd)
	return 0

def shred(infile, outfile, length):
	cmd = "shred.sh in={} out={} length={} overwrite overlap={} threads=1 2>>shred.log".format(infile, outfile, length, length-1)
	runcmd(cmd)
	return 0

def runcmd(cmd):
	#sys.stderr.write("[system] {}\n".format(cmd))
	return system(cmd)

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Shreds transcript sequences by gene to create short FASTA reads. Used as a first step in generating a mappability table for genes based on read alignments.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('fasta', type=str, help="FASTA reference of transcript sequences")
parser.add_argument('ref', type=str, help="refFlat annotation for grouping transcripts into genes.")
parser.add_argument('-p', type=int, default=1, help="Number of processes to launch for shredding.")
parser.add_argument('-l', type=int, default=75, help="Read length to generate")
parser.add_argument('-o', type=str, default=None, help="Stub for output file name. Default is to use the input fasta name stub")
parser.add_argument('-z', action="store_const", const=True, default=False, help="Gzip final shred file.")
parser.add_argument('-b', '--bbmap', action="store_const", const=True, default=False, 
	help="When finished creating shreds, map the shreds back to the reference with bbmap")
parser.add_argument('-t', type=int, default=1, help="Number of threads for bbmap")
parser.add_argument('-s', action="store_const", const=True, default=False, 
	help="Use this to bypass the shredding and go straight to alignment. Expects input to be a shredded file. Enables -b automatically. In this mode the 'ref' argument is expected to be the reference FASTA to which we will align the shreds.")
parser.add_argument('-q', '--quantify', action="store_const", const=True, default=False, 
	help="Shred, align and quantify the relationships between genes.")

parser.add_argument('--just-alignment', action="store_const", const=True, default=False, 
	help="Forces the program to jump straight to alignment. Assumes the expected shred reads file already exists.")

parser.add_argument('--just-quantification', action="store_const", const=True, default=False, 
	help="Forces the program to jump straight to quantificaiton of the alignments. Expects intended files already exist.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

