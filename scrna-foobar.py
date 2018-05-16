#!/usr/bin/python
#==============================================================================
# scrna-foobar.py
#
# Shawn Driscoll
# 20171207
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Not totally sure of all this is gonna do. Basically I think I want to be able
# to take in a bunch of per-cell demuxed single-cell reads, pool them, align
# them to a transcriptome and then generate UMI counts per gene. This code
# can alternativly filter the single-cell reads by either mapping them or by 
# kmer. 
#
# Something interesting. I've found many reads in one of our experiments that
# map perfectly to one or two genes but when aligned to the genome they 
# also map perfectly or nearly perfectly to many random locations and the aligner
# reports very low mapq - something that should probably be discarded. I hadn't
# seen this before since transcript sequence tends to be very specific. It's 
# possible, to be on the safe side, that it may be wise to filter reads
# for extreme multi-mappers against the genome:
# bowtie2 -x ~/alind/ucsc/mm10/bwt2/mm10 --local -p 12 -U <reads>.fastq | samtools view -hS -q 10 -o bowtie.sam -
# reformat.sh in=bowtie.sam out=filtered.fastq
#==============================================================================

import sys
import argparse
import math
import re
import traceback
import os.path
from collections import defaultdict
from time import localtime, time
import subprocess as sp
from os import system, unlink
from multiprocessing import cpu_count, Process, JoinableQueue, Queue, current_process, Lock
import hashlib
from random import random

# from igraph import *
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

HOME = os.path.expanduser("~")
SAMHITS = "/home/colfax/coding/python/scrna-samhits.py"
GTOT = "/home/colfax/coding/python/bam-galign-to-talign-fast.py"
GSNAP = "/home/colfax/opt/gsnap/bin/gsnap"

#==============================================================================
# main
#==============================================================================

def main(args):

	err = False

	##
	## confirm all files exist. check for target file.
	if not os.path.isfile(args.fin):
		error_message("Input file {} does not exist".format(args.fin))
		err = True
	
	##
	## input file is a text file with a list of all of the files to process
	## check that those files exist
	flist = []
	with open(args.fin, "r") as fin:
		for szl in fin:
			aln = szl.strip().split("\t")
			flist.append(aln[0])
			if not os.path.isfile(aln[0]):
				error_message("File listed in index, {}, does not exist".format(aln[0]))
				err = True

	if err:
		return 1

	# quantification mode
	doquant = (args.bbmap_quantify is not None) or (args.bowtie2_quantify is not None) or (args.bowtie_quantify is not None) or (args.genome_quantify is not None) or (args.hisat2_quantify is not None)

	if doquant:
		if args.r is None:
			error_message("You must specify -r when in quantification mode")
			err = True
		else:
			if not os.path.isfile(args.r):
				error_message("refFlat annotation you provided does not exist")
				err = True
			

	##
	## all good?
	##

	tblfile = "batch_samples.tsv"
	if not doquant:
		# remove current tblfile because this code will append to it
		if os.path.isfile(tblfile):
			unlink(tblfile)

	# loop in chunks
	bclist = []
	if args.chunk_size > 0:
		num_files = len(flist)
		steps = int(math.floor(num_files/args.chunk_size))
		if steps==0:
			steps = 1
		istart = 0
		iend = 0
		for i in range(steps):
			if i == (steps-1):
				iend = num_files
			else:
				iend = istart + args.chunk_size
			# slice the file lists
			flistHat = flist[istart:iend]
			# process the file list
			rres = core(flistHat, args, tblfile)
			
			istart = iend
			
			bclist += rres
	
	##
	# if we used STAR lets unload the reference from memory
	if (args.genome_quantify is not None):
		cmd = "star --genomeLoad Remove --genomeDir {}".format(args.genome_quantify)
		runcmd(cmd)
	
	if (args.star is not None):
		cmd = "star --genomeLoad Remove --genomeDir {}".format(args.star)
		runcmd(cmd)
	
	##
	# if we did quantification then we need to pull all of the quantification
	# files together.
	##
	if doquant:
		t0 = time()
		message("Combining quantification files")
		
		# check for quantification files
		err = False
		quantList = []
		cell_reads = {}
		for bcid in bclist:
			fname = "{}.umi_count".format(bcid)
			if not os.path.isfile(fname):
				error_message("quantification {} is missing".format(fname))
				err = True
			else:
				quantList.append(fname)
				cellid = re.sub("\.umi_count", "", fname)
				cell_reads[cellid] = 0

		if len(quantList) > 1:
			# load all of the files and export in the format of kallisto pseduo batch
			
			# get complete gene list
			glist = set()		
			for fname in quantList:
				with open(fname, "r") as fin:
					tmp = []
					for szl in fin:
						if szl[0]=="#":
							# check if this comment line contains the total read count that scrna-samhits.py includes
							r = re.search("total_reads=([0-9]+)", szl)
							if r:
								cellid = re.sub("\.umi_count", "", fname)
								cell_reads[cellid] = int(r.group(1))
							continue
						aln = szl.strip().split("\t")
						tmp.append(aln[0])
				
				glist.update(tmp)
			
			# now we have a full list. make a dict to convert gene names
			# to indices
			glist = list(glist)
			dgenes = {}
			for i in range(len(glist)):
				dgenes[glist[i]] = i
			
			# build output
			outlines = "#gene\tcell\tumi\n"
			
			for i in range(len(quantList)):
				fname = quantList[i]
				with open(fname, "r") as fin:
					for szl in fin:
						if szl[0]=="#":
							continue
						aln = szl.strip().split("\t")
						if float(aln[1]) > 0:
							# record this hit
							outlines += "\t".join(map(str, [dgenes[aln[0]], i, aln[1]]))
							outlines += "\n"
			
			# now we can write the three files
			with open("scrna.tcc", "w") as fout:
				fout.write(outlines)
			
			with open("scrna.cells", "w") as fout:
				for f in quantList:
					cellid = re.sub("\.umi_count", "", f)
					fout.write("{}\t{}\n".format(cellid, cell_reads[cellid]))
			
			with open("scrna.genes", "w") as fout:
				for gname in glist:
					fout.write("{}\n".format(gname))
			
			# remove all of the umi_count files
			for f in quantList:
				unlink(f)
		
		sys.stderr.write("{} sec\n".format(time()-t0))	

	return 0

##
# main code of the entire pipeline. it will completely process any
# files in 'flist'.
def core(flist, args, tblfile):
	
	pool = []
	taskQ = None
	tlock = None
	
	doquant = (args.bbmap_quantify is not None) or (args.bowtie_quantify is not None) or (args.bowtie2_quantify is not None) or (args.genome_quantify is not None) or (args.hisat2_quantify is not None)
	
	tmpfile = hashlib.md5(args.fin).hexdigest()
	tmpfile1 = tmpfile + ".1"
	tmpfile2 = tmpfile + ".2"
	if args.F:
		## fasta and not fastq
		tmpfile1 += ".fasta"
		tmpfile2 += ".fasta"
	else:
		tmpfile1 += ".fastq"
		tmpfile2 += ".fastq"

	if doquant:
		tmpfile2 = tmpfile + ".sam"

	if os.path.isfile(tmpfile1):
		unlink(tmpfile1)
	
	if os.path.isfile(tmpfile2):
		unlink(tmpfile2)

	message("Pooling input files")
	t0 = time()
	rres = pool_files(tmpfile1, flist, args.F, args.samplerate)
	sys.stderr.write("{} sec\n".format(time()-t0))
	
#	if args.samplerate < 1:
#		message("Subsampling reads")
#		t0 = time()
#		rres = subsample_reads(tmpfile1, args.F, args.samplerate)
#		sys.stderr.write("{} sec\n".format(time()-t0))
	
	if args.bowtie2_filter is not None:
		t0 = time()
		message("Filtering against the genome with bowtie2")
		tmpfile3 = "{}.3.fastq".format(tmpfile)
		threads = args.t
		if threads==0:
			threads = cpu_count()/2
			
		rres = bowtie2_filter(tmpfile1, tmpfile3, args.bowtie2_filter, args.bowtie2_mapq, threads)
		# move tmpfile3 to tmpfile1
		cmd = "mv {} {}".format(tmpfile3, tmpfile1)
		system(cmd)
		sys.stderr.write("{} sec\n".format(time()-t0))
	
	t0 = time()
	if args.bbduk is not None:
		# bbduk filter
		message("Filtering with bbduk")
		rres = bbduk_filter(tmpfile1, tmpfile2, args.bbduk, args.k, args.t)
	elif args.bbmap is not None:
		message("Filtering with bbmap")
		rres = bbmap_filter(tmpfile1, tmpfile2, args.bbmap, args.alignment_identity, args.t)
	elif args.star is not None:
		message("Filtering reads against the genome with STAR")
		rres = star_filter(tmpfile1, tmpfile2, args)
	elif args.bbmap_quantify is not None:
		message("Mapping reads with bbmap")
		rres = bbmap_align(tmpfile1, tmpfile2, args.bbmap_quantify, args.alignment_identity, args.t)
	elif args.bowtie_quantify is not None:
		message("Mapping reads with bowtie for quantification")
		rres = bowtie_align(tmpfile1, tmpfile2, args.bowtie_quantify, args.alignment_identity, args.t)
	elif args.bowtie2_quantify is not None:
		message("Mapping reads with bowtie2 for quantification")
		rres = bowtie2_align(tmpfile1, tmpfile2, args.bowtie2_quantify, args.alignment_identity, args.t)
	elif args.genome_quantify is not None:
		message("Aligning with STAR and translating to transcriptome hits")
		rres = star_quantify(tmpfile1, tmpfile2, args)
	elif args.hisat2_quantify is not None:
		message("Aligning with HISAT2 and translating to transcriptome hits")
		rres = hisat2_quantify(tmpfile1, tmpfile2, args)

	sys.stderr.write("{} sec\n".format(time()-t0))

	if doquant:
		##
		# deal with quantification
		message("Indexing cells within alignments")
		t0 = time()
		bc_offsets = parse_barcodes_sam(tmpfile2)
		sys.stderr.write("{} sec\n".format(time()-t0))
		
		##
		## setup processing queue
		##
		
		pool = []
		taskQ = JoinableQueue()
		tlock = Lock()
		threads = args.t
		if threads==0:
			threads = cpu_count()/2
			
		for i in range(threads):
			p = Process(target=quant_child, args=(taskQ, args, SAMHITS, tlock,))
			p.daemon = True
			p.start()
			pool.append(p)
		
		##
		## visit each barcode, write alignments out to a temp file and quantify
		## using 'scrna-samhits.py'
		##
		t0 = time()
		message("Quantifying individual cells")
		
		fin = open(tmpfile2, "r")
		
		# grab the header lines
		sam_header = ""
		for szl in fin:
			if szl[0]=="@":
				sam_header += szl
			else:
				break
		
		for bc in bc_offsets.keys():
			cellfile = "{}.{}.cell.sam".format(tmpfile, bc)
			with open(cellfile, "w") as fout:
				
				#fout.write("@HD\tVN:1.0\tSO:unsorted\n")
				fout.write(sam_header)
				
				for offset in bc_offsets[bc]:
					fin.seek(offset)
					szl = fin.readline()
					##
					# strip off the cell id so that the quantification program 
					# doesn't get confused.
					aln = szl.strip().split("\t")
					tmp = aln[0].split(":")
					aln[0] = ":".join(tmp[0:(len(tmp)-1)])
					szl = "\t".join(aln) + "\n"
					# write line to file
					fout.write(szl)
				
			# quantify
			quantfile = "{}.umi_count".format(bc)
			#fout = open(quantfile, "w")
			#cmd = "{} {} {}".format(SAMHITS, args.r, cellfile)
			#sys.stderr.write("CMD: {}\n".format(cmd))
			#p1 = sp.Popen(cmd.split(), stdout=fout)
			#p1.wait()
			#fout.close()
			#unlink(cellfile)
			
			taskQ.put([cellfile, quantfile])
		
		# close the alignments
		fin.close()
		
		# wait here for child threads to complete and wrap it all up
		for p in pool:
			taskQ.put(None)
		
		taskQ.join()
		
		for p in pool:
			p.join()
		
		sys.stderr.write("{} sec\n".format(time()-t0))
		
		unlink(tmpfile1)
		unlink(tmpfile2)
	
	else:
		##
		# now we parse tmpfile2 to identify the cells and their umis
		message("Indexing cells in filtered reads")
		t0 = time()
		bc_offsets, bc_umi = parse_barcodes_and_umis(tmpfile2, args.F)
		sys.stderr.write("{} sec\n".format(time()-t0))
	
		##
		# now we can write individual files
		message("Writing individual filtered cell reads and UMI files")
		t0 = time()
		rres = write_filtered_files(tmpfile2, args.F, bc_offsets, bc_umi, args.min_umi_freq)
		sys.stderr.write("{} sec\n".format(time()-t0))
	
		unlink(tmpfile1)
		unlink(tmpfile2)
	
		##
		# write the batch table
		szout = ""
		for cellid in rres.keys():
			num_reads = rres[cellid][3]
			
			if num_reads >= args.depth_filter:
				rres[cellid].append("PASS")
			else:
				rres[cellid].append("LOWDEPTH")
				
			szout += "\t".join(map(str, rres[cellid]))
			szout += "\n"
	
		with open(tblfile, "a") as fout:
			fout.write(szout)	
	
	return bc_offsets.keys()

##
# quant_child
# this is a worker function for the process queue which will deal with quantifying
# individual cell files in the background.
def quant_child(tasks, args, cmd_path, lock):
	# get name
	name = current_process().name
	
	# loop over tasks
	while True:
		task = tasks.get()
		if task is None:
			# done
			tasks.task_done()
			break
		
		if os.path.isfile(task[0]):
			
			# step one is to convert alignments to transcriptome based. pass the input file to a 
			# temporary file in the process
			tmpfile = hashlib.md5(task[0]).hexdigest() + ".tmpfile.sam"
			cmd = "{} --quiet --filter-mapq {} --filter-overlap {} -o {}".format(GTOT, args.filter_mapq, args.filter_overlap, tmpfile)
			if args.f_stranded:
				cmd += " --f-stranded"
			elif args.r_stranded:
				cmd += " --r-stranded"
			
			cmd += " {} {}".format(args.r, task[0])

			lock.acquire()
			sys.stderr.write("{}: {}\n".format(name, cmd))
			# release lock
			lock.release()
					
			system(cmd)
		
			# run quantification. task is a two item list with the 
			# input file name and target output file name
			#cmd = "{} --quiet {} {} > {}".format(cmd_path, args.r, task[0], task[1])
			cmd = "{} --quiet {} {} > {}".format(cmd_path, args.r, tmpfile, task[1])
			# get lock
			lock.acquire()
			sys.stderr.write("{}: {}\n".format(name, cmd))
			# release lock
			lock.release()
			system(cmd)

			# remove input file
			unlink(task[0])
			unlink(tmpfile)
			
		else:
			lock.acquire()
			sys.stderr.write("{}: input file does not exist ({})\n".format(name, task[0]))
			lock.release()
			
		tasks.task_done()
	
	return 0

def gzip_child(tasks, lock):
	
	while True:
		task = tasks.get()
		if task is None:
			tasks.task_done()
			break
		
		# task is the file name
		if os.path.isfile(task):
			cmd = "gzip {}".format(task)
			system(cmd)

		tasks.task_done()
	
	return 0
	

def write_filtered_files(infile, fasta, bco, bcu, min_freq):
	##
	# output file names will be the barcode id plus ".filtered.[fastq|fasta]"
	
	##
	# process for compressing read files
	taskQ = JoinableQueue()
	tlock = Lock()
	threads = 2
	pool = []
	for i in range(threads):
		p = Process(target=gzip_child, args=(taskQ, tlock,))
		p.daemon = True
		p.start()
			
	dout = {}
	fin = open(infile, "r")
	
	for cellid in bco.keys():
		progress_message("Exporting cell {} ({} reads)".format(cellid, len(bco[cellid])))
		offsets = bco[cellid]
		szout = ""
		szumi = ""
		dout[cellid] = [cellid, "", "", 0, 0]
		
		for o in offsets:
			fin.seek(o)
			
			szl = fin.readline()
			# parse barcodes

			bc = parse_barcodes(szl.strip())
			umi = bc[0]
			if bcu[cellid][umi] >= min_freq:
				# we can print this one
				dout[cellid][3] += 1
				
				szumi += "{}\n".format(umi)
				
				# we can keep this read
				# strip off the cell id
				aln = szl.strip().split(":")				
				szout += ":".join(aln[0:(len(aln)-1)])
				szout += "\n"
				
				if fasta:
					while True:
						# read from file until we hit another read
						szl = fin.readline()
						if szl[0]==">":
							break
							
						szout += szl
				else:
					idx = 3
					while idx > 0:
						idx -= 1
						szout += fin.readline()
			else:
				# count umis that are dropped by filter
				dout[cellid][4] += 1
				
		outfile = "{}.filtered.".format(cellid)
		if fasta:
			outfile += "fasta"
		else:
			outfile += "fastq"
		
		dout[cellid][2] = outfile
		with open(outfile, "w") as fout:
			fout.write(szout)
		
		# compress it in the background
		taskQ.put(outfile)
		
		outfile = "{}.filtered.umi".format(cellid)
		dout[cellid][1] = outfile
		with open(outfile, "w") as fout:
			fout.write(szumi)
	
	fin.close()
		
	sys.stderr.write("\n")

	sys.stderr.write("waiting for gzip child process to complete...\n")
	for i in range(threads):
		taskQ.put(None)

	taskQ.join()
	
	for p in pool:
		p.join()

	return dout
	
def parse_barcodes_sam(infile):
	
	# parse the sam file, containing alignments from multiple cells, 
	# to index where all the reads from individual cells are located.
	
	bcindex = defaultdict(list)
	offset = 0
	
	with open(infile, "r") as fin:
		for szl in fin:
			if szl[0]=="@":
				offset += len(szl)
				continue
			
			aln = szl.strip().split("\t")
			bc = parse_barcodes(aln[0])
			bcindex[bc[1]].append(offset)
			
			offset += len(szl)
	
	return bcindex
	

def parse_barcodes_and_umis(infile, fasta):
	
	# dict for barcodes. this will hold the line offset position for all reads
	# that belong to a barcode
	dbc = {}
	# dict for barcodes. this will hold a dict per barcode key'd by umis and
	# holding their counts
	dbcu = {}
	
	offset = 0
	idx = 0
	
	with open(infile, "r") as fin:
		
		for szl in fin:
			idx += 1
			if (szl[0]==">" and fasta) or ((not fasta) and idx == 1):
				# read name line
				bc = parse_barcodes(szl.strip())
				cellid = bc[1]
				umi = bc[0]

				if cellid not in dbc:
					dbc[cellid] = []
					dbcu[cellid] = {}

				if umi not in dbcu[cellid]:
					dbcu[cellid][umi] = 0
				
				dbc[cellid].append(offset)
				dbcu[cellid][umi] += 1
	
			if idx > 3:
				idx = 0
			
			offset += len(szl)
	
	return dbc, dbcu

##
# this function will open all files and pool them into a single file. it's going to 
# add the file's name to the read names so that we can demux the reads again
# after filtering.
def pool_files(outfile, flist, fasta, samplerate):
	
	num_reads = 0
	write_read = False
	
	if os.path.isfile(outfile):
		unlink(outfile)
	
	if fasta:
		tmpfile = "{}.tmp.fasta".format(outfile)
	else:
		tmpfile = "{}.tmp.fastq".format(outfile)
	
	if os.path.isfile(tmpfile):
		unlink(tmpfile)
	
	if samplerate < 1:
		message("Subsampling rate: {}".format(samplerate))
	
	for f in flist:
		progress_message("pooling {} ".format(f))
		# pick apart the file name
		tmp = f.split("/")
		fname = tmp[-1]
		tmp2 = fname.split(".")
		fstub = tmp2[0]
		
		# gzipped?
		if re.search("\.gz$", f):
			p1 = sp.Popen("gunzip -c {}".format(f).split(), stdout=sp.PIPE)
			fin = p1.stdout
		else:
			fin = open(f, "r")
		
		with open(tmpfile, "a") as fout:
			idx = 0
			write_read = False
			for szl in fin:
				
				idx += 1
				if (fasta and szl[0]==">") or ((not fasta) and idx==1):
					rx = random()
					
					if rx < samplerate:
						write_read = True
						num_reads += 1
						# read name line. append the stub for this file so we can track it
						# drop any whitespace delim stuff too
						aln = szl.strip().split()
						szl = aln[0] + ":{}\n".format(fstub)
					else:
						write_read = False
				
				if write_read:
					fout.write(szl)
					
				if idx > 3:
					idx = 0
	
		fin.close()
	
		
	sys.stderr.write("\n")
	sys.stderr.write("Parsed {} total reads\n".format(num_reads))
	
	##
	## shuffle the reads
	##
	if False:
		sys.stderr.write("Shuffling...\n")
		cmd = "shuffle.sh in={} out={}".format(tmpfile, outfile)
		system(cmd)
		unlink(tmpfile)
	
	cmd = "mv {} {}".format(tmpfile, outfile)
	runcmd(cmd)
	
	return num_reads

def get_filetype(fname):

	tmp = str.lower(fname).split(".")
	n = len(tmp)

	gz = tmp[-1] == "gz"

	if gz:
		ext = tmp[-2]
	else:
		ext = tmp[-1]

	if ext == "fasta" or ext == "fa":
		fmt = "fasta"
	elif ext == "fastq" or ext == "fq":
		fmt = "fastq"

	return [fmt, gz]

def subsample_reads(infile, fasta, rate):
	if fasta:
		outfile = "{}.sub.fasta".format(infile)
	else:
		outfile = "{}.sub.fastq".format(infile)
	
	cmd = "reformat.sh in={} out={} samplerate={}".format(infile, outfile, rate)
	rres = runcmd(cmd)
	if os.path.isfile(outfile):
		# overwrite input with subsampled reads
		cmd = "mv {} {}".format(outfile, infile)
		runcmd(cmd)
	
	return 0

##
# this function returns the cell barcode and umi barcode from a read name string
def parse_barcodes(sz):
	tmp = sz.split(":")
	return [tmp[-2], tmp[-1]]

##
# runs bbduk to filter reads against a reference keeping any read with 
# at least a single kmer match
def bbduk_filter(infile, outfile, ref, k, threads):

	cmd = "bbduk.sh in={} ref={} k={} hdist=0 mm=f overwrite=t threads={}".format(infile, ref, k, threads)
	cmd += " outm={}".format(outfile)
	
		
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	
	return 0

##
# runs bbmap to filter reads keeping only those that pass some minimum 
# alignment identity
def bbmap_filter(infile, outfile, ref, identity, threads):
	
	minid = identity
	bbmap_opts = "maxindel=1 minid={} local=t overwrite=t usejni=f threads={}".format(identity, threads)
	cmd = "bbmap.sh in={} path={} idfilter={} outm={} {}".format(infile, ref, identity, outfile, bbmap_opts)
		
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	
	return 0


##
# runs bbmap to filter reads keeping only those that pass some minimum 
# alignment identity
def bbmap_align(infile, outfile, ref, identity, threads):
	
	minid = identity
	bbmap_opts = "maxindel=1 minid={} local=t overwrite=t usejni=f ambig=all threads={}".format(identity, threads)
	cmd = "bbmap.sh in={} path={} idfilter={} outm={} {}".format(infile, ref, identity, outfile, bbmap_opts)
		
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	
	return 0

def bowtie2_filter(infile, outfile, ref, mapq, threads):
	
	# run bowtie2
	cmd = "bowtie2 -x {} -U {} -p {} --local | samtools view -F 0x4 -hS -q {} - | reformat.sh in=stdin.sam out={}".format(ref, infile, threads, mapq, outfile)
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	
	# convert back to fastq
	#cmd = "reformat.sh in={}.sam out={}".format(outfile, outfile)
	#sys.stderr.write("CMD: {}\n".format(cmd))
	#system(cmd)
	
	# remove the sam file
	#samfile = "{}.sam".format(outfile)
	#if os.path.isfile(samfile):
	#	unlink(samfile)
	
	return 0

def bowtie2_align(infile, outfile, ref, identity, threads):	
	identity = 1-identity
	cmd = "bowtie2 -x {} -U {} -p {} -k 200 --mp 1,1 --np 1 --score-min L,0,-{} > {}".format(ref, infile, threads, identity, outfile)
	return runcmd(cmd)

def bowtie_align(infile, outfile, ref, identity, threads):
	
	eval = 70 # bowtie default
	mean_qual = 22
	
	# check read length
	with open(infile, "r") as fin:
		szl = fin.readline()
		# read the query string
		szl = fin.readline()
		rl = len(szl.strip())
		message("Detected read length {}".format(rl))
		nm = math.floor(rl * (1-identity))
		eval = nm * mean_qual
		message("Using -e {}".format(eval))
		
	
	cmd = "bowtie -p {} -e {} -k 500 --best --chunkmbs 128 -S {} {} > {}".format(threads, int(math.floor(eval)), ref, infile, outfile)
	return runcmd(cmd)
	

def star_quantify(infile, outfile, args):
	
	mm = 1 - args.alignment_identity
	threads = args.t
	
	if threads==0:
		threads = cpu_count()/2
		
	tmpfile = "Aligned.out.sam"
	cmd = "star-aln -t {} -n {} --no-unal --sam --load-and-keep --sensitive --un -k 1000".format(threads, mm)
	cmd += " {} {}".format(args.genome_quantify, infile)
	rres = runcmd(cmd)
	
	if False:
		# check to see if we have a bbmap index in this same location
		args.genome_quantify = re.sub("\/$", "", args.genome_quantify)
		tmp = (args.genome_quantify).split("/")
		tmp[len(tmp)-1] = "bbmap"
		tmp = "/".join(tmp)
		if os.path.isdir(tmp):
			# ok, attemp to map the unaligned reads with bbmap
			runcmd("mv Unmapped.out.mate1 Unmapped.fq")
			cmd = "bbrun bbmap -in Unmapped.fq -path {} -ambig all -slow -secondary -maxsites 1000 -maxindel 500000 -intronlen 21 -sam 1.3 -threads {} -outm Aligned.bbmap.sam".format(tmp, threads)
			runcmd(cmd)
			runcmd("grep -v '^@' Aligned.bbmap.sam >> Aligned.out.sam")
			
			if os.path.isfile("Aligned.bbmap.sam"):
				unlink("Aligned.bbmap.sam")
	
	kills = ["Log.out", "Log.progress.out", "Log.final.out", "Unmapped.out.mate1", "Unmapped.fq", "SJ.out.tab"]
	for fname in kills:
		if os.path.isfile(fname):
			unlink(fname)
	
	if os.path.isfile(tmpfile):
		# ok good deal. now we need to send this alignment to the genome to transcriptome tool
		cmd = "{} --filter-mapq {} --filter-overlap {} -o {}".format(GTOT, args.filter_mapq, args.filter_overlap, outfile)
		if args.f_stranded:
			cmd += " --f-stranded"
		elif args.r_stranded:
			cmd += " --r-stranded"
		
		cmd += " {} {}".format(args.r, tmpfile)
				
		rres = runcmd(cmd)
		
		unlink(tmpfile)
	
	return 0


def hisat2_quantify(infile, outfile, args):
	
	mm = 1 - args.alignment_identity
	threads = args.t
	
	if threads==0:
		threads = cpu_count()/2
	
	# see if there is an 'ss.txt' file in the index folder
	index_parts = os.path.split(args.hisat2_quantify)
	ss_file = index_parts[0] + "/ss.txt"
	
	tmpfile = "Aligned.out.sam"
	cmd = "hisat2 -p {} -x {} -k 200 --score-min L,0.0,-{} -U {} ".format(threads, args.hisat2_quantify, mm, infile)
	
	if os.path.isfile(ss_file):
		cmd += "--known-splicesite-infile {} ".format(ss_file)
	
	#cmd += "> {}".format(tmpfile)
	cmd += "> {}".format(outfile)
	
	rres = runcmd(cmd)

	##
	## 20180227 seeing if passing this off to the multi-threaded part of the pipeline is 
	## any more efficient.
	if False:
		
		if os.path.isfile(tmpfile):
			# ok good deal. now we need to send this alignment to the genome to transcriptome tool
			cmd = "{} --filter-mapq {} --filter-overlap {} -o {}".format(GTOT, args.filter_mapq, args.filter_overlap, outfile)
			if args.f_stranded:
				cmd += " --f-stranded"
			elif args.r_stranded:
				cmd += " --r-stranded"
			
			cmd += " {} {}".format(args.r, tmpfile)
					
			rres = runcmd(cmd)
			
			unlink(tmpfile)
	
	return 0


##
## added 1/5/2018
## I updated my installation of gsnap and found that it does a significantly better job
## at finding alignments than STAR in certain cases and with mostly default settings.
def gsnap_quantify(infile, outfile, args):
	
	mm = 1 - args.alignment_identity
	threads = args.t
	
	if mm > 0.1:
		sys.stderr.write("Warning: gsnap doesn't like mismatch ratios greater than 0.1. Chaning your setting to 0.1\n")
		mm = 0.1
	
	if threads==0:
		threads = cpu_count()/2
	
	# parse path to gsnap index. make sure to drop any trailing '/' so that the following 
	# split command doesn't mess up
	if re.search("\/$", args.genome_quantify):
		tmp = re.sub("\/$", "", args.genome_quantify)
		args.genome_quantify = tmp
	
	# split path up to extract the genome name from the end
	index_tmp = (args.genome_quantify).split("/")
	gsnap_genome = index_tmp[-1]
	gsnap_dir = "/".join(index_tmp[0:(len(index_tmp)-1)])
	
	# make temporary alignment file name
	tmpfile = infile + ".gsnap.sam"
	cmd = GSNAP + " -N 1 -t {} -m {} -A sam -o {}".format(threads, mm, tmpfile)
	cmd += " -D {} -d {}".format(gsnap_dir, gsnap_genome)
	if re.search("\.gz$", infile):
		cmd += " --gunzip"
	
	cmd += " {}".format(infile)
#	cmd = "star-aln -t {} -n {} --no-unal --sam --load-and-keep -k 200".format(threads, mm)
#	cmd += " {} {}".format(args.genome_quantify, infile)

	rres = runcmd(cmd)
	if os.path.isfile(tmpfile):
		# ok good deal. now we need to send this alignment to the genome to transcriptome tool
		cmd = "{} --filter-mapq {} --filter-overlap {} -o {}".format(GTOT, args.filter_mapq, args.filter_overlap, outfile)
		if args.f_stranded:
			cmd += " --f-stranded"
		elif args.r_stranded:
			cmd += " --r-stranded"
		
		cmd += " {} {}".format(args.r, tmpfile)
		
		rres = runcmd(cmd)
		
		unlink(tmpfile)
	
	return 0


def star_filter(infile, outfile, args):

	mm = 1 - args.alignment_identity
	threads = args.t
	
	if threads==0:
		threads = cpu_count()/2
		
	tmpfile = "Aligned.out.sam"

	cmd = "star-aln -t {} -n {} --no-unal --sam --load-and-keep -k 10".format(threads, mm)
	cmd += " {} {}".format(args.star, infile)
	cmd += " seedSearchStartLmax=15 winAnchorMultimapNmax=200"
	rres = runcmd(cmd)
	
	if os.path.isfile(tmpfile):
		# convert to fastq
		cmd = "samtools view -S -q {} -F 0x100 {} | reformat.sh in=stdin.sam out={}".format(args.filter_mapq, tmpfile, outfile)
		runcmd(cmd)
		unlink(tmpfile)
	
	return 0
	

def runcmd(cmd):
	sys.stderr.write("CMD: {}\n".format(cmd))
	system(cmd)
	return 0

##
# filters reads based on umi frequency. all distinct umis are counted and then 
# reads are written to 'outfile' that have UMIs which occur at least the 
# specified minimum number of times.
def umi_filter(infile, outfile, min_count, fmt, gz):

	# if data is gzipped then first we will undo that
	if gz:
		cmd = "gunzip {}".format(infile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)
		infile = re.sub("\.gz", "", infile)

	if fmt == "fasta":
		rres = umi_filter_fasta(infile, outfile, min_count)
	else:
		rres = umi_filter_fastq(infile, outfile, min_count)

	if gz:
		cmd = "gzip {}".format(infile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)
		cmd = "gzip {}".format(outfile)
		sys.stderr.write("CMD: {}\n".format(cmd))
		system(cmd)

	return rres


def umi_filter_fasta(infile, outfile, min_count):

	idx = 0
	dumi = defaultdict(int)
	keep_read = False
	keep_count = 0
	drop_count = 0

	##
	# first count the UMIs
	with open(infile, "r") as fin:

		for szl in fin:
			if szl[0] == ">":
				# read name, get barcodes
				bc = parse_barcodes(szl.strip())
				dumi[bc[1]] += 1

	##
	# read back through and print out the reads passing the filter. also write
	# the umi file.
	szout = ""
	umiout = ""
	with open(infile, "r") as fin:

		for szl in fin:

			if szl[0] == ">":
				keep_read = False
				bc = parse_barcodes(szl.strip())
				if dumi[bc[1]] >= min_count:
					umiout += "{}\n".format(bc[1])
					keep_read = True
					keep_count += 1
				else:
					drop_count += 1

			if keep_read:
				szout += szl

	with open(outfile, "w") as fout:
		fout.write(szout)
	
	tmp = outfile.split(".")
	tmp[-1] = "umi"
	umifile = ".".join(tmp)
	with open(umifile, "w") as fout:
		fout.write(umiout)
	

	return [keep_count, drop_count]


def umi_filter_fastq(infile, outfile, min_count):

	idx = 0
	dumi = defaultdict(int)
	keep_read = False
	keep_count = 0
	drop_count = 0

	##
	# first count the UMIs
	with open(infile, "r") as fin:

		for szl in fin:
			idx += 1
			if idx == 1:
				# read name, get barcodes
				bc = parse_barcodes(szl.strip())
				dumi[bc[1]] += 1

			if idx == 4:
				idx = 0

	##
	# read back through and print out the reads passing the filter
	szout = ""
	umiout = ""
	idx = 0
	with open(infile, "r") as fin:

		for szl in fin:
			idx += 1
			if idx == 1:
				keep_read = False
				bc = parse_barcodes(szl.strip())
				if dumi[bc[1]] >= min_count:
					umiout += "{}\n".format(bc[1])
					keep_read = True
					keep_count += 1
				else:
					drop_count += 1

			if keep_read:
				szout += szl

			if idx == 4:
				idx = 0

	with open(outfile, "w") as fout:
		fout.write(szout)

	tmp = outfile.split(".")
	tmp[-1] = "umi"
	umifile = ".".join(tmp)
	with open(umifile, "w") as fout:
		fout.write(umiout)

	return [keep_count, drop_count]


def progress_message(sz, last=False):
	tmp = [" " for i in range(80)]
	sys.stderr.write("\r{}".format("".join(tmp)))
	sys.stderr.write("\r> {}".format(sz))
	if last:
		sys.stderr.write("\n")
	return 0


def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def error_message(sz):
	sys.stderr.write("Error: {}\n".format(time_string(), sz))

def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def print_exception():
	exc_type, exc_value, exc_traceback = sys.exc_info()
	print "*** print_tb:"
	traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
	print "*** print_exception:"
	traceback.print_exception(exc_type, exc_value, exc_traceback,
	                          limit=2, file=sys.stdout)
	print "*** print_exc:"
	traceback.print_exc()
	print "*** format_exc, first and last line:"
	formatted_lines = traceback.format_exc().splitlines()
	print formatted_lines[0]
	print formatted_lines[-1]
	print "*** format_exception:"
	print repr(traceback.format_exception(exc_type, exc_value,
	                                      exc_traceback))
	print "*** extract_tb:"
	print repr(traceback.extract_tb(exc_traceback))
	print "*** format_tb:"
	print repr(traceback.format_tb(exc_traceback))
	print "*** tb_lineno:", exc_traceback.tb_lineno	

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('fin', type=str, help="Text file index of read files to process, one file with full path per line.")

parser.add_argument("-F", action="store_const", const=True, default=False, 
	help="Input is FASTA format and not FASTQ")

parser.add_argument("-t", action="store", default=0, type=int, 
	help="Threads for aligner and BBTools programs. If set to zero it will use all_cores/2.")

parser.add_argument('--chunk-size', action="store", default=100, type=int, 
	help="Number of cell files to process per chunk. If set to 0 then all are processed at once.")

parser.add_argument("--bowtie2-filter", action="store", default=None, type=str, 
	help="Filter reads against the genome using bowtie2. Useful for filtering out excessive multi-mapping reads to the genome.")

parser.add_argument("--bowtie2-mapq", action="store", default=20, type=int, 
	help="MAPQ filter when using --bowtie2-filter")

myg = parser.add_mutually_exclusive_group(required=True)
myg.add_argument("--bbduk", type=str, default=None, 
	help="FASTA reference to filter reads against with bbduk.")
myg.add_argument("--bbmap", type=str, default=None, 
	help="Path to bbmap index to filter reads against")
myg.add_argument("--star", type=str, default=None,
	help="Path to star genome index to filter reads against")
myg.add_argument("--bbmap-quantify", type=str, default=None, 
	help="Path to bbmap transcriptome index to quantify UMI against. You must also provide an annotation with -r")
myg.add_argument("--bowtie2-quantify", type=str, default=None, 
	help="Path to bowtie2 transcriptome index to quantify UMI against. You must also provide an annotation with -r")
myg.add_argument("--bowtie-quantify", type=str, default=None, 
	help="Path to bowtie transcriptome index to quantify UMI against. You must also provide an annotation with -r")
myg.add_argument("--genome-quantify", type=str, default=None, 
	help="Path to STAR genome index. Reads are aligned to genome and translated to transcriptome hits and finally quantified against refFlat annotation specified with -r")
myg.add_argument("--hisat2-quantify", type=str, default=None, 
	help="Path to HISAT2 genome index. Reads are aligned to genome and translated to transcriptome hits and finally quantified against refFlat annotation specified with -r")


parser.add_argument("-r", type=str, default=None, 
	help="refFlat annotation to use if using --bbmap-quantify")

parser.add_argument("--samplerate", type=float, default=1.0, 
	help="Subsample reads per cell at this rate. This is executed in the bbmap or bbduk call.")

parser.add_argument("-k", type=int, default=31, 
	help="Kmer size for read filtering against a reference.")

parser.add_argument("--alignment-identity", type=float, default=0.8, 
	help="Alignment identity as ratio of total read length")
	
parser.add_argument("--depth-filter", type=int, default=10000, 
	help="Minimum read depth filter. Applied after other filtering steps. This adds a flag in the 'batch_samples.txt' file.")

parser.add_argument('--min-umi-freq', type=int, default=1, 
	help="Minimum frequency of a UMI within a cell for output to file.")

parser.add_argument('--filter-mapq', type=int, default=4, 
	help="Used with --genome-quantify. Minimum MAPQ for genome alignments.")

parser.add_argument('--filter-overlap', type=float, default=0.75, 
	help="Used with --genome-quantify. Minimum overlap between alignment and transcript. If >= 1 then it is interpreted as a number of bases. If 0 < x < 1 then it is interpreted as a ratio of the read length")

strand_group = parser.add_mutually_exclusive_group(required=False)
strand_group.add_argument('--f-stranded', action="store_const", const=True, default=False, 
	help="Used with --genome-quantify. Keep only alignments that have the same strand orientation as the target transcript")
strand_group.add_argument('--r-stranded', action="store_const", const=True, default=False, 
	help="Used with --genome-quantify. Keep only alignments that have the opposite strand orientation as the target transcript")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		print_exception()
		sys.exit(1)

