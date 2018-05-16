#!/usr/bin/python
#==============================================================================
# species-split-ab.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Split RNA-Seq reads between two species via genome alignment and 
# score based sdorting. 
#==============================================================================

import sys, argparse, math, re, os
from os.path import isfile, expanduser
from hashlib import md5
import subprocess as sp
import multiprocessing as mp
import timeit
from time import localtime, time
import traceback
from Basics import  utils
import gzip

# import numpy as np

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")

CHUNK_SIZE = 1000000
NUM_PHASE_ONE_PROC = 2
NUM_PHASE_TWO_PROC = 2
NUM_PHASE_THREE_PROC = 1

# sam flags
SAM_PAIRED = 0x1
SAM_PROPERLY_PAIRED = 0x2
SAM_UNAL = 0x4
SAM_MATE_UNAL = 0x8
SAM_REVERSED = 0x10
SAM_MATE_REVERSED = 0x20
SAM_FIRST_MATE = 0x40
SAM_SECOND_MATE = 0x80
SAM_SECONDARY = 0x100

# sam alignment fields
SAM_QNAME = 0
SAM_FLAG = 1
SAM_RNAME = 2
SAM_POS = 3
SAM_MAPQ = 4
SAM_CIGAR = 5
SAM_RNEXT = 6
SAM_PNEXT = 7
SAM_TLEN = 8
SAM_SEQ = 9
SAM_QUAL = 10
SAM_FIRST_ATTR = 11

# python for converting sam format to fastq
BAM2FASTQ = "{}/coding/python/bam2fastq.py".format(HOME)
BAM2FASTQ += " -o {} {}"

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	
	line_count = 0
	read_count = 0
	chunk_count = 0
	buffer_handle = None
	buffer_fname = ""
	buffer_count = 0
	tmpname = ""
	p = None
	phase1_pool = []
	phase2_pool = []
	phase3_pool = []
	phase1_tasks = mp.JoinableQueue()
	phase2_tasks = mp.JoinableQueue()
	phase3_tasks = mp.JoinableQueue()
	results = mp.Queue()
	
	lock = mp.Lock()
	
	#
	# main loop is going to parse the input reads file and write the reads
	# back out to a buffer file. if nothing else were to happen the main loop
	# would essentially just split the input reads file into N files where
	# N = ceiling(total_reads/CHUNK_SIZE). when a chunk is completed we will 
	# pass it into a queue where the initial set of alignments happen.
	#
	
	#
	# check args
	#
	
	if not isfile(args.reads):
		error_message("input reads cannot be found {}".format(args.reads))
		return 1
	
	# check indexes
	
	#
	# initalize stuff
	#
	
	tmpname = md5(args.reads).hexdigest()
	buffer_fname = "{}_{}.fq".format(tmpname, buffer_count)
	buffer_handle = open(buffer_fname, "w")
	run_t0 = timeit.default_timer()
	
	#
	# load STAR indexes into shared memory
	#

	t0 = time()
	message("Loading STAR indexes into shared memory to speed things up a bit")
	indexes = [args.indexA, args.indexB]
	for idx in indexes:

		cmd = "star --genomeLoad Remove --genomeDir {}".format(idx)
		sys.stderr.write("CMD: {}\n".format(cmd))
		p = sp.Popen(cmd.split())
		p.wait()

		cmd = "star --genomeLoad LoadAndExit --genomeDir {}".format(idx)
		sys.stderr.write("CMD: {}\n".format(cmd))
		p = sp.Popen(cmd.split())
		p.wait()
	
	sys.stderr.write("{} seconds\n".format(time()-t0))

	# 
	# start a child process that will deal with the individual read file chunks
	#
	
	message("Starting processes")
	for i in range(NUM_PHASE_ONE_PROC):
		p = mp.Process(target=phase_one_worker, args=(phase1_tasks, phase2_tasks, args, lock))
		p.daemon = True
		p.start()
		phase1_pool.append(p)
	
	for i in range(NUM_PHASE_TWO_PROC):
		p = mp.Process(target=phase_two_worker, args=(phase2_tasks, phase3_tasks, args, lock))
		p.daemon = True
		p.start()
		phase2_pool.append(p)

	for i in range(NUM_PHASE_THREE_PROC):
		p = mp.Process(target=phase_three_worker, args=(phase3_tasks, args, lock))
		p.daemon = True
		p.start()
		phase3_pool.append(p)
	
	# 
	# get the reads open and get going
	#
	message("Opening reads")
	fin = None
	if re.search("\.gz$", args.reads):
		# reads are gzipped
#		message("opening gzipped reads")
#		cmd = "gunzip -c {}".format(args.reads)
#		p1 = sp.Popen(cmd.split(), stdout=sp.PIPE)
#		fin = p1.stdout
		fin = gzip.open(args.reads, "r")
	else:
		fin = open(args.reads, "r")
	
	if fin is None:
		error_message("No input stream was established")
		return 1
	
	for szl in fin:
		line_count += 1

		# write current line to buffer
		buffer_handle.write(szl)

		if line_count==4:
			# count read
			read_count += 1
			chunk_count += 1
			line_count = 0
		
		if chunk_count >= CHUNK_SIZE:
			message("parsed {} reads".format(read_count))
			
			# close the buffer and start a new one
			buffer_handle.close()
			#
			# pass to queue
			phase1_tasks.put(buffer_fname)
			# make a new one
			buffer_count += 1
			buffer_fname = "{}_{}.fq".format(tmpname, buffer_count)
			buffer_handle = open(buffer_fname, "w")
			chunk_count = 0
			
	
	# close input
	fin.close()

	# close the buffer and pass to tasks
	buffer_handle.close()
	phase1_tasks.put(buffer_fname)
	
	message("finished splitting reads. waiting for children...")
	
	# insert stop
	for p in phase1_pool:
		phase1_tasks.put(None)
	# wait for phase one to finish before putting in stops
	phase1_tasks.join()

	# send stops to phase 2
	for p in phase2_pool:
		phase2_tasks.put(None)

	phase2_tasks.join()

	# send stops to phase 2
	for p in phase3_pool:
		phase3_tasks.put(None)
		
	phase3_tasks.join()

	for p in phase1_pool:
		p.join()
	for p in phase2_pool:
		p.join()
	for p in phase3_pool:
		p.join()
	
	# loop through results queue and build final files

	# we need some output files. from the queue we will have a bunch of stuff to 
	# pool
	labs = (args.l).split(",")

	a_unique = "{}_{}_unique.fq.gz".format(args.o, labs[0])
	a_ambig = "{}_{}_ambig.sam".format(args.o, labs[0])
	b_unique = "{}_{}_unique.sam".format(args.o, labs[1])
	b_ambig = "{}_{}_ambig.sam".format(args.o, labs[1])
	ab_ambig = "{}_ambig.sam".format(args.o)
	low_score = "{}_lowscore.sam".format(args.o)

	a_ambigR = "{}_{}_ambig.fq.gz".format(args.o, labs[0])
	b_uniqueR = "{}_{}_unique.fq.gz".format(args.o, labs[1])
	b_ambigR = "{}_{}_ambig.fq.gz".format(args.o, labs[1])
	ab_ambigR = "{}_ambig.fq.gz".format(args.o)
	low_scoreR = "{}_lowscore.fq.gz".format(args.o)

	ab_unal = "{}_unal.fq".format(args.o)
	b_unal = "{}_{}_unal.fq".format(args.o, labs[1])

	#
	# now with all of these things pooled we can align the B unaligned
	# reads to A and generate the final A unique and AB unaligned 
	# results
	#
	
	# 
	# we also have to translate b_unique, b_ambig, ab_ambig, and a_ambig 
	# from sam to FASTQ
	#
	message("Running final alignments to {}".format(labs[0]))
	rres = phase_three(args, b_unal, ab_unal, a_unique, low_score)
	
	# 
	# some of these files need to be translated to fastq
	message("Translating SAM files to FASTQ and compressing")
	flist = [a_ambig, b_unique, b_ambig, ab_ambig, low_score]
	flistR = [a_ambigR, b_uniqueR, b_ambigR, ab_ambigR, low_scoreR]
	for i in range(len(flist)):
		try:
			os.unlink(flistR[i])
		except:
#			print_exception()
			pass
			
		#cmd = "reformat.sh in={} out={} overwrite 1>/dev/null 2>/dev/null".format(flist[i], flistR[i])
		cmd = BAM2FASTQ.format(flistR[i], flist[i])
		rres = utils.runcmd(cmd, args.v)
		if rres[0] == 0:
			os.unlink(flist[i])
		else:
			error_message("reformat failed on {}".format(flist[i]))
	
	# need to compress ab_unal
	try:
		os.unlink("{}.gz".format(ab_unal))
	except:
		pass
	cmd = "gzip {}".format(ab_unal)
	utils.runcmd(cmd, args.v)
	
	# count to make stats
	message("counting sort results and writing to stats file")
	flist = [a_unique, a_ambigR, b_uniqueR, b_ambigR, ab_ambigR, low_scoreR, "{}.gz".format(ab_unal)]
	cmd = "{}/coding/perl/bb-count.pl {} >{}.tsv".format(HOME, " ".join(flist), tmpname)
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		error_message("bb-count failed!")
	
	# open the file up and insert percentages
	fin = open("{}.tsv".format(tmpname), "r")
	aln = []
	olines = []
	tot_reads = 0
	for szl in fin:
		aln = szl.strip().split("\t")
		tmp = re.sub(args.o + "_", "", aln[0])
		lout = [args.o, re.sub("\.fq\.gz$", "", tmp), aln[1]]
		tot_reads += int(aln[1])
		olines.append(lout)
	fin.close()
	
	fout = open("{}.tsv".format(args.o), "w")
	# write header
	fout.write("sample\tgroup\tnum_reads\tpercent_reads\ttotal_reads\n")
	for lout in olines:
		lout += ["{:0.2f}".format(int(lout[-1])*100.0/tot_reads), tot_reads]
		fout.write("\t".join(map(str, lout)) + "\n")
	
	fout.close()
	
	try:
		os.unlink("{}.tsv".format(tmpname))
	except:
#		print_exception()		
		pass		
	
	#
	# UNload STAR indexes 
	#

	if not args.no_unload:
		message("Unloading STAR indexes from shared memory")
		cmd = "star --genomeLoad Remove --genomeDir {}".format(args.indexA)
		utils.runcmd(cmd, verbose=True)
		cmd = "star --genomeLoad Remove --genomeDir {}".format(args.indexB)
		utils.runcmd(cmd, verbose=True)
	
	run_time = timeit.default_timer()-run_t0
	run_min = run_time*1.0/60
	run_sec = math.ceil((run_min-math.floor(run_min))*60)
	run_min = math.floor(run_min)
	message("Completed in {:d} minutes and {:d} seconds".format(int(run_min), int(run_sec)))
	
	for f in ["Aligned.out.sam", "Log.out", "Log.progress.out", "_STARtmp"]:
		try:
			os.unlink(f)
		except:
			pass
			
	
	return 0


#
# this is the 'main' for when we're just aligning the sort results from a 
# previous run
def main2(args):
	
	labs = (args.l).split(",")
	
	#
	# first figure out if the expected files are present
	#
	stub = args.o
	# expected files names
	a_unique = "{}_{}_unique.fq.gz".format(stub, labs[0])
	a_ambig = "{}_{}_ambig.fq.gz".format(stub, labs[0])
	b_unique = "{}_{}_unique.fq.gz".format(stub, labs[1])
	b_ambig = "{}_{}_ambig.fq.gz".format(stub, labs[1])
	flist = [a_unique, a_ambig, b_unique, b_ambig]
	
	a_bams = ["{}_{}_unique.bam".format(stub, labs[0]), "{}_{}_ambig.bam".format(stub, labs[0])]
	b_bams = ["{}_{}_unique.bam".format(stub, labs[1]), "{}_{}_ambig.bam".format(stub, labs[1])]
	a_final_bam = "{}_{}.bam".format(stub, labs[0])
	b_final_bam = "{}_{}.bam".format(stub, labs[1])
	
	# check files
	missing_files = False
	for fname in flist:
		if not isfile(fname):
			missing_files = True
			error_message("Missing expected sort file {}".format(fname))
		else:
			message("Found {}".format(fname))
	
	if missing_files:
		error_message("Some files were missing. Bailing out")
		return 1
		
	#--------------------------------------------------------------------------
	#
	# ALIGNMENT TO SPECIES A
	#
	#--------------------------------------------------------------------------
	cmd0 = "{}/coding/perl/star-aln.pl -t {} -n {} --load-and-keep --cord-sort".format(HOME, args.p, args.final_map_n)
	if args.bowtie2_A != "":
		cmd0 += " --bowtie-ref {} --bowtie-local".format(args.bowtie2_A)
	
	cmd = cmd0 + " --rgid unique {} {} 2>/dev/null".format(args.indexA, a_unique)
	utils.runcmd(cmd, True)
	
	# when finished we need to rename the alignments
	if isfile("Aligned.sortedByCoord.out.bam"):
		cmd = "mv Aligned.sortedByCoord.out.bam {}".format(a_bams[0])
		utils.runcmd(cmd, True)
	else:
		error_message("Expected output from STAR is not here!")
		return 1
	
	# map the ambig file
	cmd = cmd0 + " --rgid ambig {} {} 2>/dev/null".format(args.indexA, a_ambig)
	utils.runcmd(cmd, True)
	
	utils.runcmd("star --genomeLoad Remove --genomeDir {}".format(args.indexA), True)

	if not isfile("Aligned.sortedByCoord.out.bam"):
		error_message("Expected output from STAR is not here!")
		return 1		
	
	# merge these
	cmd = "samtools merge -f -@ {} {} {} {}".format(args.p, a_final_bam, a_bams[0], "Aligned.sortedByCoord.out.bam")
	utils.runcmd(cmd, True)

	#--------------------------------------------------------------------------
	#
	# ALIGNMENT TO SPECIES B
	#
	#--------------------------------------------------------------------------
	cmd0 = "{}/coding/perl/star-aln.pl -t {} -n {} --load-and-keep --cord-sort".format(HOME, args.p, args.final_map_n)
	if args.bowtie2_B != "":
		cmd0 += " --bowtie-ref {} --bowtie-local".format(args.bowtie2_B)
	
	cmd = cmd0 + " --rgid unique {} {} 2>/dev/null".format(args.indexB, b_unique)
	utils.runcmd(cmd, True)
	
	# when finished we need to rename the alignments
	if isfile("Aligned.sortedByCoord.out.bam"):
		cmd = "mv Aligned.sortedByCoord.out.bam {}".format(b_bams[0])
		utils.runcmd(cmd, True)
	else:
		error_message("Expected output from STAR is not here!")
		return 1
	
	# map the ambig file
	cmd = cmd0 + " --rgid ambig {} {} 2>/dev/null".format(args.indexB, b_ambig)
	utils.runcmd(cmd, True)
	
	utils.runcmd("star --genomeLoad Remove --genomeDir {}".format(args.indexB), True)
	
	if not isfile("Aligned.sortedByCoord.out.bam"):
		error_message("Expected output from STAR is not here!")
		return 1		
	
	# merge these
	cmd = "samtools merge -f -@ {} {} {} {}".format(args.p, b_final_bam, b_bams[0], "Aligned.sortedByCoord.out.bam")
	utils.runcmd(cmd, True)
	

	#--------------------------------------------------------------------------
	#
	# CLEAN UP
	#
	#--------------------------------------------------------------------------
	flist = ["Aligned.sortedByCoord.out.bam", "SJ.out.tab", "Log.final.out", "Log.progress.out", "Log.out", a_bams[0], b_bams[0]]
	for fname in flist:
		try:
			message("removing {}".format(fname))
			os.unlink(fname)
		except:
			print_exception()
			pass
	
	# all done
	
	return 0
	

#==============================================================================
# defs
#==============================================================================

def dump_file_to_handle(fname, handle):
	try:
		fin = open(fname, "r")
		for szl in fin:
			handle.write(szl)
		fin.close()
	except:
		print_exception()
		pass
	
	return


def phase_one_worker(task_queue, results_queue, args, lock):
	name = mp.current_process().name
	message_mp("Phase one process started", name, lock)
	rres = None
	labs = (args.l).split(",")
	stub = ""
	
	while True:
		task = task_queue.get()
		if task is None:
			task_queue.task_done()
			break
		
		# do stuff
		message_mp("P1: starting {}".format(task), name, lock)
		stub = re.sub("\.fq$", "", task)
#		try:
#			rres = phase_one(args, task, stub, name, lock)
#			# remove the reads file
#			os.unlink(task)
#		except Exception, e:
#			message_mp("something went wrong", name, lock)
#			print_exception()
#			continue
		
		rres = phase_one(args, task, stub, name, lock)
		os.unlink(task)
		
		if rres[0] is None:
			# something went wrong
			message_mp(rres[1], name, lock)
			# don't put anything into the result queue
			task_queue.task_done()
		else:
			rres.append(stub)
			message_mp("P1: finished {}".format(task), name, lock)
			results_queue.put(rres)
			task_queue.task_done()
	
	message_mp("P1: shutting down", name, lock)	
	
	return
		

def phase_two_worker(task_queue, results_queue, args, lock):
	name = mp.current_process().name
	message_mp("Phase two process started", name, lock)
	labs = (args.l).split(",")
	rres = None
	
	while True:
		task = task_queue.get()
		if task is None:
			task_queue.task_done()
			break
		
		# we should get a list of file names back
		#message_mp("received {}".format(",".join(task)), name, lock)
		# 0: aligned to B
		# 1: reads of 0 as aligned to A
		# 2: unaligned reads from B alignment (pass through)
		# 3: stub for files
		
		message_mp("P2: received {}".format(task[3]), name, lock)
		rres = phase_two(task[0], task[1], task[3], args)
		if rres[0] is None:
			# something went wrong
			message_mp(rres[1])
			task_queue.task_done()
		else:
			# pass through the unaligned reads from phase 1
			rres.append(task[2])
			results_queue.put(rres)
			task_queue.task_done()
			
		# drop files needed only for this phase 
		os.unlink(task[0])
		os.unlink(task[1])
	
		message_mp("P2: finished {}".format(task[3]), name, lock)
	
	message_mp("P2: shutting down", name, lock)
	
	return

#
# this worker receives the split data file names from phase2 and concats the 
# results into sample level files
def phase_three_worker(task_queue, args, lock):
	name = mp.current_process().name
	message_mp("Phase three process started", name, lock)
	labs = (args.l).split(",")
	
	# file names
	b_unique = "{}_{}_unique.sam".format(args.o, labs[1])
	b_ambig = "{}_{}_ambig.sam".format(args.o, labs[1])
	a_ambig = "{}_{}_ambig.sam".format(args.o, labs[0])
	ab_ambig = "{}_ambig.sam".format(args.o)
	low_score = "{}_lowscore.sam".format(args.o)
	b_unal = "{}_{}_unal.fq".format(args.o, labs[1])
	stub = ""

	# file list in order of the items that will show up from tasks	
	flist = [b_unique, b_ambig, a_ambig, ab_ambig, low_score, b_unal]
	handles = [open(flist[i], "w") for i in range(len(flist))]
	
	while True:
		item = task_queue.get()
		if item is None:
			task_queue.task_done()
			break
		
		# derive stub from the first file name
		r = re.search("^([^\_]+\_[0-9]+)", item[0])
		if r:
			stub = r.group(1)
		else:
			stub = item[0]
		
		message_mp("P3: received {}".format(stub), name, lock)
		
		try:	
			# dump contents of the files from phase 2 into the cumulative files	
			i = 0
			while i < len(item):
				dump_file_to_handle(item[i], handles[i])
				# delete the item file
				os.unlink(item[i])					
				i += 1
				
		except:
			print_exception()			
			pass

		message_mp("P3: finished {}".format(stub), name, lock)
		
		task_queue.task_done()
	
	# close the handles
	for i in range(len(handles)):
		handles[i].close()		
		
	message_mp("P3: shutting down", name, lock)
	
	return		

def phase_one(args, reads, stub, proc_name, lock):

	use_bowtie2 = False
	labs = (args.l).split(",")
	pid = None
	childs = []
	num_proc = max([1, args.p/NUM_PHASE_ONE_PROC])
	
	#
	# star settings
	star_settings = {
		'genomeLoad': "LoadAndKeep",
		'runMode': "alignReads", 
		'runThreadN': num_proc,
		'outFilterScoreMinOverLread': 0.6, 
		'outFilterMatchNminOverLread': 0, 
		'alignEndsType': "Local", 
		'outSAMtype': "SAM",
		'outFilterMultimapScoreRange': 1, 
		'outFilterMultimapNmax': 1000, 
		'outFilterMismatchNmax': 1000, 
		'outFilterMismatchNoverLmax': 1, 
		'outFilterMismatchNoverReadLmax': 1, 
		'outFilterScoreMin': 0, 
		'outFilterScoreMinOverLread': 0.6, 
		'outFilterMatchNmin': 0, 
		'outFilterMatchNminOverLread': 0.5, 
		'outQSconversionAdd': -31 if args.q else 0,
		'outFileNamePrefix': "{}_".format(stub), 
		'outSAMattributes': "NH HI AS nM NM MD"
	}
	
	# file name list

	b_aligned = "{}_{}_aln.sam".format(stub, labs[1])
	b_aligned_reads = "{}_{}_aln.fq".format(stub, labs[1])
	b_unal_reads = "{}_{}_unal.fq".format(stub, labs[1])
	ba_aligned = "{}_{}_{}_aln.sam".format(stub, labs[1], labs[0])

	# star files
	star_files = ["Aligned.out.sam", "Unmapped.out.mate1", "Log.out", "SJ.out.tab", "Log.final.out", "Log.progress.out"]
	for i in range(len(star_files)):
		star_files[i] = "{}_{}".format(stub, star_files[i])

	#'readFilesIn': reads, 
	#'genomeDir': index,
	#'outSAMunmapped': "None" if unal_out else "Within", 
	#'outReadsUnmapped': "Fastx" if unal_out else "None",
	#'outSAMattrRGline': "ID:{}".format(rgid), 

	if args.bowtie2_B != "" and args.bowtie2_A != "":
		use_bowtie2 = True

	# align to B with star. we need unaligned reads to be sent to a fastq file
	# and the aligned need to be sent to both SAM and FASTQ
	message_mp("mapping to {} with STAR".format(labs[1]), proc_name, lock)
	cmd = "star --outSAMattrRGline ID:{} ".format(labs[1])
	cmd += " --readFilesIn {} --genomeDir {}".format(reads, args.indexB)
	cmd += " --outReadsUnmapped Fastx --outSAMunmapped None"
	cmd += " {} 1>/dev/null".format(collapse_star_settings(star_settings))
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		# none
		return [None, "Failed to align to B with STAR"]
	
	
	# unmapped reads are in 'Unmapped.out.mate1'
	# alignments are in 'Aligned.out.sam'
	
	if use_bowtie2:
		# align  unaligned reads to B with bowtie2. send unmapped reads out to a fastq file.
		# these reads are those that fail to map to B and must be mapped to A later.
		message_mp("mapping STAR unaligned reads with bowtie2", proc_name, lock)
		cmd = "bowtie2 -x {} -U {}_Unmapped.out.mate1 --un {} ".format(args.bowtie2_B, stub, b_unal_reads)
		cmd += " -p {} --local --rg-id {}".format(num_proc, labs[1])
		cmd += " 2>/dev/null | samtools view -SF 0x4 - >>{}_Aligned.out.sam".format(stub)
		rres = utils.runcmd(cmd, args.v)
		if rres[0] != 0:
			return [None, "Failed to align to B with bowtie2"]
			
	else:
		# no bowtie2 so 'Unmapped.out.mate1' will be the unaligned set that fail to align to B
		# and should 
		cmd = "mv {}_Unmapped.out.mate1 {}".format(stub, b_unal_reads)
		rres = utils.runcmd(cmd, args.v)
		if rres[0] != 0:
			return [None, "Failed to rename Unmapped.out.mate1 after aligning to B with STAR"]
		
		
	#
	# need to do two things which have be executed in series 
	# 1. filter alignments for only primary alignments
	# 2. copy those alignments to a fastq file
	#

	cmd = "samtools view -hSF 0x100 -o {} {}_Aligned.out.sam".format(b_aligned, stub)
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		return [None, "Failed to filter out secondary alignments from B alignents"]

	#cmd = "reformat.sh in={} out={} overwrite 1>/dev/null 2>/dev/null".format(b_aligned, b_aligned_reads)
	cmd = BAM2FASTQ.format(b_aligned_reads, b_aligned)
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		return [None, "Failed to translate B alignments to FASTQ"]

		
	#
	# b_aligned_reads are now mapped to A. if we're using bowtie2 then the unaligned
	# reads have to be exported. otherwise we keep them in the SAM file for the 
	# sorting stage
	#

	# remove previous star files
	for f in star_files:
		try:
			os.unlink(f)
		except:
			pass

	message_mp("mapping {} aligned reads to {} with STAR".format(labs[1], labs[0]), proc_name, lock)

	cmd = "star --outSAMattrRGline ID:{} ".format(labs[0])
	cmd += " --readFilesIn {} --genomeDir {}".format(b_aligned_reads, args.indexA)
	
	if use_bowtie2:
		cmd += " --outReadsUnmapped Fastx --outSAMunmapped None"
	else:
		cmd += " --outReadsUnmapped None --outSAMunmapped Within"
	
	# if the base qualities were phred64 they have been modified now so we can disable this
	star_settings['outQSconversionAdd'] = 0	
	cmd += " {} 1>/dev/null".format(collapse_star_settings(star_settings))
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		return [None, "Failed to align B aligned reads to A with STAR"]
	
	if use_bowtie2:
		# align the unmapped reads to A
		message_mp("mapping STAR unaligned reads with bowtie2", proc_name, lock)
		cmd = "bowtie2 -x {} -U {}_Unmapped.out.mate1 -p {} --local --rg-id {} --no-head 2>/dev/null >>{}_Aligned.out.sam".format(args.bowtie2_A, stub, num_proc, labs[0], stub)
		rres = utils.runcmd(cmd, args.v)
		if rres[0] != 0:
			return [None, "Failed to align B aligned reads to A with bowtie2"]

	# alignments have to be filtered for primary alignments only. 
	cmd = "samtools view -hSF 0x100 -o {} {}_Aligned.out.sam".format(ba_aligned, stub)
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		return [None, "Failed to filter out secondary alignments from B_A aligned reads"]
	
	# now we want to pool b_aligned and ba_aligned, sort by read name and go through the 
	# B/A sorting process. so do i want to push the pair of file names out to some 
	# child processing queue so we can move forward? is that too much since this 
	# is an alignment thread?
	
	star_files += [b_aligned_reads]
	for f in star_files:
		try:
			os.unlink(f)
		except:
			pass
	
	return [b_aligned, ba_aligned, b_unal_reads]

#
# in this phase we pool the Baligned and BAaligned reads, sort by 
# read name and then evaluate B/A assignment
def phase_two(Baligned, BAaligned, stub, args):

	labs = (args.l).split(",")
	dlabs = {}
	for i in range(len(labs)):
		dlabs[labs[i]] = i
		
	# output files for sorted reads
	b_unique = "{}_{}_unique.sam".format(stub, labs[1])
	b_ambig = "{}_{}_ambig.sam".format(stub, labs[1])
	a_ambig = "{}_{}_ambig.sam".format(stub, labs[0])
	ambig = "{}_ambig.sam".format(stub)
	low_score = "{}_lowscore.sam".format(stub)
	pool = "{}_pool.sam".format(stub)
	pool_sorted = "{}_sorted.sam".format(stub)
	
	# minimum score gap for assignment. within two mismatch is ambiguous
	min_score_gap = args.mp*2
	
	# for parsing the sam lines
	rname = ""
	rname_last = ""
	ab_aln = [[], []]
	
	# pool the alignment files and then sort by read name
	try:
		fout = open(pool, "w")
		fin = open(Baligned, "r")
		for szl in fin:
			if szl[0]=="@":
				continue
			fout.write(szl)
		fin.close()
		fin = open(BAaligned, "r")
		for szl in fin:
			if szl[0]=="@":
				continue
			fout.write(szl)
		fin.close()
		fout.close()
	except:
		print_exception()		
		pass
	
	# sort
	cmd = "sort -k1,1 {} > {}".format(pool, pool_sorted)
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		return [None, "Failed to sort pooled B and B_A alignments"]
	
#	try:
#		os.unlink(Baligned)
#		os.unlink(BAaligned)
#	except:
#		pass
	
	# we have to open that pooled, sorted file
	try:
		
		fin = open(pool_sorted, "r")
		# open all of the output files
		fout_b_unique = open(b_unique, "w")
		fout_b_ambig = open(b_ambig, "w")
		fout_a_ambig = open(a_ambig, "w")
		fout_ambig = open(ambig, "w")
		fout_low_score = open(low_score, "w")
		
		for szl in fin:
			# these are sam alignments
			szl = szl.strip()
			aln = szl.split("\t")
			
			rname = aln[0]
			
			if rname != rname_last and rname_last != "":
				if len(ab_aln[0]) > 0 or len(ab_aln[1]) > 0:
					# we have reads to process!
					selAB = read_select_ab(ab_aln, args)					
					
					if selAB is not None:
						# we have a sort. if it was None then no alignments went in
						# check score gap
						if selAB[2] is None:
							# unique assignment
							if selAB[0]==0:
								warning_message("for some reason a B aligned read uniquely sorted to A")
							elif selAB[0]==1:
								# assign to B
								fout_b_unique.write(ab_aln[1][0] + "\n")
						
						else:
							if selAB[2] >= min_score_gap:
								# assign to ambiguous for A or B
								if selAB[0]==0:
									fout_a_ambig.write(ab_aln[0][0] + "\n")
								else:
									fout_b_ambig.write(ab_aln[1][0] + "\n")
							else:
								# gap is too small so this is ambig
								fout_ambig.write(ab_aln[0][0] + "\n")								
					
					# clear buffer
					ab_aln = [[], []]
			
			# find RG tag value and put the read into the appropriate bin
			if (int(aln[1]) & 0x4)==0:
				#
				# read is aligned, read group tag
				#
				r = re.search("RG\:Z\:([A-Z0-9a-z]+)", szl)
				if r:
					index = dlabs[r.group(1)]
					# append
					ab_aln[index].append(szl)
				else:
					warning_message("unable to find RG tag in\nREAD: {}".format(szl))
			
			rname_last = rname
			
		# 
		# handle final bundle of reads
		#
		if len(ab_aln[0]) > 0 or len(ab_aln[1]) > 0:
			# we have reads to process!
			selAB = read_select_ab(ab_aln, args)
			
			if selAB is not None:
				# we have a sort. if it was None then no alignments went in
				
				# check threshold pass flag
				if selAB[3]==1:
					# passed threshold, continue...
				
					# check score gap
					if selAB[2] is None:
						# unique assignment
						if selAB[0]==0:
							warning_message("for some reason a B aligned read uniquely sorted to A")
						elif selAB[0]==1:
							# assign to B
							fout_b_unique.write(ab_aln[1][0] + "\n")
					
					else:
						if selAB[2] >= min_score_gap:
							# assign to ambiguous for A or B
							if selAB[0]==0:
								fout_a_ambig.write(ab_aln[0][0] + "\n")
							else:
								fout_b_ambig.write(ab_aln[1][0] + "\n")
						else:
							# gap is too small so this is ambig
							fout_ambig.write(ab_aln[0][0] + "\n")								
				
				else:
					# did not pass score threshold so this has to be passed to unaligned
					fout_low_score.write(ab_aln[0][0] + "\n")
					
				
		# all done
		
		fin.close()
		fout_b_unique.close()
		fout_b_ambig.close()
		fout_a_ambig.close()
		fout_ambig.close()
		fout_low_score.close()
		
		os.unlink(pool)
		os.unlink(pool_sorted)
		
	except Exception, e:
		print_exception()
		pass
	
	# this list gets bounded to phase_three
	return [b_unique, b_ambig, a_ambig, ambig, low_score]


# 
# arguments are the args from command line parsing and three file names.
# we'll align b_uanl to A and then concatenate unaligned reads to 
# ab_unal and send aligned to a_unique
def phase_three(args, b_unal, ab_unal, a_unique, low_score):
	
	use_bowtie2 = False
	labs = (args.l).split(",")
	
	#
	# star settings
	star_settings = {
		'genomeLoad': "LoadAndKeep",
		'runMode': "alignReads", 
		'runThreadN': args.p,
		'outFilterScoreMinOverLread': 0.6, 
		'outFilterMatchNminOverLread': 0, 
		'genomeDir': args.indexA,
		'readFilesIn': b_unal,
		'outSAMunmapped': "None",
		'outReadsUnmapped': "Fastx",
		'alignEndsType': "Local", 
		'outSAMtype': "SAM",
		'outFilterMultimapScoreRange': 1, 
		'outFilterMultimapNmax': 1000, 
		'outFilterMismatchNmax': 1000, 
		'outFilterMismatchNoverLmax': 0.1, 
		'outFilterMismatchNoverReadLmax': 1, 
		'outFilterScoreMin': 0, 
		'outFilterScoreMinOverLread': 0.6, 
		'outFilterMatchNmin': 0, 
		'outFilterMatchNminOverLread': 0.5, 
		'outQSconversionAdd': 0
	}

	if args.bowtie2_B != "" and args.bowtie2_A != "":
		use_bowtie2 = True

	# align to B with star. we need unaligned reads to be sent to a fastq file
	# and the aligned need to be sent to both SAM and FASTQ
	
	cmd = "star " + collapse_star_settings(star_settings) + " 1>/dev/null"
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		# none
		return [None, "Failed to align final reads to A with STAR"]
	
	# unmapped reads are in 'Unmapped.out.mate1'
	# alignments are in 'Aligned.out.sam'
	
	
	if use_bowtie2:
		# align  unaligned reads to A with bowtie2. 

		# temp file for unaligned reads if using bowtie2
		a_unal = "{}_{}_unal.fq".format(args.o, labs[0])

		cmd = "bowtie2 -x {} -U Unmapped.out.mate1 --un {} ".format(args.bowtie2_A, a_unal)
		cmd += " -p {} --local --rg-id {} 2>/dev/null".format(args.p, labs[1])
		cmd += " | samtools view -SF 0x4 - >>Aligned.out.sam"
		rres = utils.runcmd(cmd, args.v)
		if rres[0] != 0:
			return [None, "Failed to align final reads to A with bowtie2"]
		
		# append unalinged
		cmd = "cat {} >> {}".format(a_unal, ab_unal)
		utils.runcmd(cmd, args.v)
		os.unlink(a_unal)
			
	else:
		# no bowtie2 so 'Unmapped.out.mate1' will be the fianl unaligned set
		cmd = "cat Unmapped.out.mate1 >> {}".format(ab_unal)
		rres = utils.runcmd(cmd, args.v)
		if rres[0] != 0:
			return [None, "Failed to rename Unmapped.out.mate1 after aligning to B with STAR"]
		
	
	# filter on alignment score 
	message("Filtering {} alignments on score".format(labs[0]))
	with open("Aligned.out.sam", "r") as fin, open("Aligned.filtered.sam", "w") as fout_good, open(low_score, "a") as fout_ls:
		for szl in fin:
			if szl[0]=="@":
				continue
			
			aln = szl.strip().split("\t")
			if (int(aln[SAM_FLAG]) & SAM_SECONDARY) != 0:
				continue
			
			aln_score = sam_score_alignment(aln, args.io, args.ie, args.mp, args.ma)
			if aln_score[1] < args.st:
				fout_ls.write(szl)
			else:
				fout_good.write(szl)
	
	# translate to FASTQ
	#cmd = "reformat.sh in=Aligned.filtered.sam out={} overwrite=t 1>/dev/null 2>/dev/null".format(a_unique)
	cmd = BAM2FASTQ.format(a_unique, "Aligned.filtered.sam")
	rres = utils.runcmd(cmd, args.v)
	if rres[0] != 0:
		return [None, "Failed to filter alignments to A and create final fastq file"]
	
	# clean up
	flist = ["Unmapped.out.mate1", "Aligned.out.sam", "Aligned.filtered.sam", "Log.progress.out", "Log.final.out", "Log.out", "SJ.out.tab", b_unal]
	for f in flist:
		try:
			os.unlink(f)
		except:
			pass
	
	
	return 0

def max_with_index(v):
	n = len(v)
	r0 = min(v)
	i0 = 0
	
	for i in range(n):
		if v[i] > r0:
			r0 = v[i]
			i0 = i
	
	return [r0, i0]

#
# this function determines which sort, A or B, an aligned read is assigned to. 
# returns the sort: A=0, B=1, the score for the winning species as well as the
# difference. the fourth item is a flag telling if the accepted final score
# is above the args.st threshold for minimum score as ratio of the aligned read's
# maximum score
def read_select_ab(buffer, args):
	
	nA = len(buffer[0])
	nB = len(buffer[1])
	
	scoresA = []
	scoresAr = []
	scoresB = []
	scoresBr = []
	
	rres = None
	
	if nA==0 and nB==0:
		return rres
	
	elif nA > 0 and nB > 0:
		# make decision
		for i in range(nA):
			tmp = sam_score_alignment(buffer[0][i].split("\t"), args.io, args.ie, args.mp, args.ma)
			scoresA.append(tmp[0])
			scoresAr.append(tmp[1])
		for i in range(nB):
			tmp = sam_score_alignment(buffer[1][i].split("\t"), args.io, args.ie, args.mp, args.ma)
			scoresB.append(tmp[0])
			scoresBr.append(tmp[1])
		
		tmp = max_with_index(scoresA)		
		bestA = tmp[0]
		bestAi = tmp[1]
		
		tmp = max_with_index(scoresB)
		bestB = tmp[0]
		bestBi = tmp[1]
		
		BAdiff = bestB-bestA
		
		if BAdiff==0:
			# same score
			rres = [None, bestB, 0, -1]
		elif BAdiff > 0:
			
			if scoresBr[bestBi] >= args.st:
				rres = [1, bestB, BAdiff, 1]
			else:
				rres = [1, bestB, BAdiff, 0]
				
		else:
			# diff must be negative
			if scoresAr[bestAi] >= args.st:
				rres = [0, bestA, abs(BAdiff), 1]
			else:
				rres = [0, bestA, abs(BAdiff), 0]			
		
	else:
		if nA > 0:
			# only alignments to A (this never would happen)
			for i in range(nA):
				tmp = sam_score_alignment(buffer[0][i].split("\t"), args.io, args.ie, args.mp, args.ma)
				scoresA.append(tmp[0])
				scoresAr.append(tmp[1])
			
			tmp = max_with_index(scoresA)
			bestA = tmp[0]
			bestAi = tmp[1]
			
			if scoresAr[bestAi] >= args.st:
				rres = [0, bestA, None, 1]
			else:
				rres = [0, bestA, None, 0]
			
		elif nB > 0:
			# only alignments to B

			for i in range(nB):
				tmp = sam_score_alignment(buffer[1][i].split("\t"), args.io, args.ie, args.mp, args.ma)
				scoresB.append(tmp[0])
				scoresBr.append(tmp[1])
			
			tmp = max_with_index(scoresB)
			bestB = tmp[0]
			bestBi = tmp[1]
				
			if scoresBr[bestBi] >= args.st:
				rres = [1, bestB, None, 1]
			else:
				rres = [1, bestB, None, 0]
	
	return rres			
	


# expecting that passed object is the list that was split from a single
# sam alignment line
def sam_aligned_length(laln):
	cigar = laln[SAM_CIGAR]
	r = re.findall("([0-9]+)[MD]", cigar)
	return sum(map(int, r))


#
# score all alignments manually using a bowtie2 local alignment score strategy.
# we have gap open/extend penalties, mismatch penalties and a  match bonus.
def sam_score_alignment(aln, gap_open, gap_ext, mismatch, match_bonus):
	
	num_match = 0
	num_mismatch = 0
	penalty = 0
	
	num_indel = 0
	num_indel_ext = 0
	
	rlen = len(aln[SAM_SEQ])
	alnlen = sam_aligned_length(aln)
	max_score = match_bonus*alnlen
	
	# get cigar
	cigar = aln[SAM_CIGAR]
	
	# ops and lens
	ops = re.findall("[0-9]+([MIDNSHP\=X])", cigar)
	lens = re.findall("([0-9]+)[MIDNSHP\=X]", cigar)
	
	# determine cigar spec
	if re.search("[\=X]", cigar):
		# sam 1.4 spec
		for i in range(len(ops)):
			if ops[i]=="=":
				num_match += int(lens[i])
			elif ops[i]=="X":
				num_mismatch += int(lens[i])
			elif ops[i]=="I" or ops[i]=="D":
				num_indel += 1
				num_indel_ext += int(lens[i])
	
	else:
		# same 1.3 spec. count matches and indels
		for i in range(len(ops)):
			if ops[i]=="M":
				num_match += int(lens[i])
			elif ops[i]=="I" or ops[i]=="D":
				num_indel += 1
				num_indel_ext += int(lens[i])
		
		# parse out mismatches from 'NM' tag
		r = re.search("NM\:i\:([0-9]+)", " ".join(aln), re.IGNORECASE)
		if r:
			num_mismatch = int(r.group(1))
		else:
			warning_message("Failed to find 'NM' tag for:\n\t{}".format(" ".join(aln)))
		
		# adjust match count for number of mismatches
		num_match -= num_mismatch
		
	# now we know mismatch count and indel number and number of extensions
	penalty = gap_open*num_indel + gap_ext*num_indel_ext + mismatch*num_mismatch
	bonus = num_match*match_bonus
	
	score = bonus-penalty
	return [score, score*1.0/max_score]
	

# 
# key value pairs
def collapse_star_settings(d):
	lkv = []
	for k in sorted(d.keys()):
		if (type(d[k]) is int) or (type(d[k] is float)):
			lkv.append("--{} {}".format(k, d[k]))
		else:
			lkv.append("--{} '{}'".format(k, d[k]))
			
	szout = " ".join(lkv)
	return szout


def time_string():
	tt = localtime()
	sz = "{}/{}/{} {:02d}:{:02d}:{:02d}".format(tt.tm_mon, tt.tm_mday, tt.tm_year, tt.tm_hour, tt.tm_min, tt.tm_sec)
	return sz

def error_message(sz):
	sys.stderr.write("[{}] Error: {}\n".format(time_string(), sz))

def warning_message(sz):
	sys.stderr.write("[{}] Warning: {}\n".format(time_string(), sz))

def message(sz):
	sys.stderr.write("[{}] {}\n".format(time_string(), sz))

def sys_message(sz):
	sys.stderr.write("[system] {}\n".format(sz))

def message_mp(sz, name, lock):
	lock.acquire()
	sys.stderr.write("[{}] {}\n".format(name, sz))
	lock.release()
	return


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

parser = argparse.ArgumentParser(description="Splits FASTQ reads into reads most likely to belong to one of two species.")

parser.add_argument('indexA', type=str, help="STAR index for species A")
parser.add_argument('indexB', type=str, help="STAR index for species B")
parser.add_argument('reads', type=str, help="Reads in FASTQ format. May be gzipped.")
parser.add_argument('--no-unload', action="store_const", const=True, default=False, 
	help="Do not unload the star genomes from shared memory")
parser.add_argument('-o', type=str, default="splitAB", help="Stub file name for all outputs. [splitAB]")
parser.add_argument('-l', type=str, default="A,B", help="Labels for the two species [A,B]")
parser.add_argument('-p', type=int, default=1, 
	help="Number of processors for alignment steps [1]")
parser.add_argument('-q', action="store_const", const=True, default=False, 
	help="Input reads have phred64 qualities [phred33 default]")
parser.add_argument('--no-q', action="store_const", const=True, default=False, 
	help="Causes the program to ignore the setting of '-q' in automated situations")
parser.add_argument('-A', '--bowtie2-A', type=str, default="", 
	help="Bowtie2 index for species A [None]")
parser.add_argument('-B', '--bowtie2-B', type=str, default="", 
	help="Bowtie2 index for species B [None]")

parser.add_argument('-v', action="store_const", const=True, default=False, 
	help="Echo all system commands as we go")

# scoring settings
group1 = parser.add_argument_group("scoring", "Alignment scoring")
group1.add_argument('--ma', type=int, default=2, help="match bonus [2]")
group1.add_argument('--mp', type=int, default=6, help="mismatch penalty [6]")
group1.add_argument('--io', type=int, default=5, help="indel open penalty [5]")
group1.add_argument('--ie', type=int, default=3, help="indel extension penalty [3]")
group1.add_argument('--st', type=float, default=0.9, 
	help="score threshold for acceptance as alignment (applied after determining genome assignment) as a ratio of a perfect score [0.9]")

# mapping only mode
group2 = parser.add_argument_group("finalmap", "Post-sort Mapping Options")
group2.add_argument('--final-mapping', action="store_const", const=True, default=False, 
	help="If the expected final output files are already present then with this option the program maps those sorted reads to their genomes.")
group2.add_argument('--final-kill-fastq', action="store_const", const=True, default=False, 
	help="After mapping the sorted reads, remove the FASTQ files (since all reads will be in the bams)")
group2.add_argument('--final-map-n', type=float, default=0.04, 
	help="All up to <int> mismatches or <num> ratio of read length mismatches [0.04]")


args = parser.parse_args()
print args
if __name__ == "__main__":

	try:
		if args.final_mapping:
			sys.exit(main2(args))
		else:
			sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

