#!/usr/bin/python
#
# bwa-quant.py
#
# porting the perl version over to python...just because
#
#

import os, sys, re
from os.path import isfile, expanduser
import argparse
from hashlib import md5
import subprocess as sp
from collections import defaultdict
import multiprocessing as mp
from math import sqrt

#==============================================================================
# Globals
#==============================================================================

HOME = expanduser("~")
BWA = HOME + "/opt/bwa-0.6.2/bwa"

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

# XA alignment fields
XA_RNAME = 0
XA_POS = 1
XA_STRAND = 4


#==============================================================================
# main
#==============================================================================

def main(args):

	pe = False
	files_exist = True
	ftmp = ""
	tasks = mp.JoinableQueue()
	p = None
	pool = []
	lock = mp.Lock()
	has_gtf = False
	files_exist = True
	sai = []
	aln_buffer = []
	frag_len_mean = 0
	frag_len_mean0 = 200
	frag_len_stdev0 = 20
	frag_len_stdev = 0
	frag_len_n = 1
	
	#
	# counters
	#
	num_parsed = 0
	num_reads = 0
	num_aligned = 0
	num_assigned = 0
	num_assigned_unique = 0
	num_assigned_multi = 0
	
	# dict for holding hits to transcripts
	dhits = defaultdict(float)
	annot = {}

	# -- -- -- check arguments
		
	if not isfile(args.index):
		# index doesn't exist
		error_message("BWA index does not exist: {}".format(args.index))
		return 1
	
	if args.g != "":
		if not isfile(args.g):
			# gtf file provided does not exist
			error_message("GTF file {} does not exist!".format(args.g))
			return 1
		# set flag for use later in the code
		has_gtf = True
	
	if len(args.reads) > 1:
		# paired end is assumed if we have two files
		pe = True
		warning_message("Detected more than one input reads file so assuming PE")
	
	# check the reads files that have been provided
	for f in args.reads:
		if not isfile(f):
			error_message("Input file {} is missing".format(f))
			files_exist = False
	
	if not files_exist:
		# at least one file didn't exist. bail out
		return 1

	
	#
	# I guess we're all good, now we may continue.
	#	
	
	#
	# create filename stub. this will be used for all 'sai' alignments as well as 
	# the 'sam' file generated via 'sampe/samse'
	ftmp = md5(args.reads[0]).hexdigest()

	#
	# fire up the child process for alignments
	#
	p = mp.Process(target=bwa_do_alignment, args=(tasks, args, ftmp, lock))
	p.daemon = True
	p.start()
	#
	# put all of the read files into a the task queue for the child process to 
	# deal with.
	for i in range(len(args.reads)):
		tasks.put([args.reads[i], i+1])
		# expected sai file result
		sai.append("{}.{}.sai".format(ftmp, i+1))
	
	# put in the stop task
	tasks.put(None)
	
	#
	# if we need to load a GTF up we can do that now while those alignments are running.
	if has_gtf:
		message("Loading {}".format(args.g))
		annot, gid2tid, gname2tid = gtf_load(args.g)
		message("Finished loading annotation. Loaded {} transcripts".format(len(annot.keys())))
		# initalize a hits table
		for tid in annot:
			dhits[tid] = 0
	else:
		#
		# no annotation but we can still grab the feature lengths from the alignments so 
		# that we can still do an 'effective length' calculation
		annot = defaultdict(int)
		gid2tid = None
		gname2tid = None
	
	#
	# wait here for alignments to complete
	message("Waiting for alignments to complete")
	tasks.join()	
	p.join()
	p.terminate()
	message("Alignments complete. Now running '{}'".format("sampe" if pe else "samse"))
	
	#
	# build the 'sampe' or 'samse' command string
	if pe:
		cmd = bwa_do_sampe(args, ftmp, sai)
	else:
		cmd = bwa_do_samse(args, ftmp, sai)
	
	# program crashes if the fifo already exists so we have to kill it if is
	# does. for some reasong 'isfile' does not return True when looking 
	# for an existing fifo so i put this into a 'try' instead.
	try:
		os.unlink("{}.sam".format(ftmp))
	except:
		# haha
		message("haha")
	
	# make fifo for alignments. we'll open this up and read it in the main 
	# loop
	os.mkfifo("{}.sam".format(ftmp))
	
	# fork out the 'samse/pe' call
	childs = []
	pid = os.fork()
	if pid==0:
		# child
		sys.stderr.write("[main|fork] {}\n".format(cmd))
		os.system(cmd)
		os._exit(0)
	else:
		childs.append(pid)
	
	# read alignments 
	message("Parsing alignments")
	rlast = rname = ""
	fin = open("{}.sam".format(ftmp), "r")
	for szl in fin:
		if szl[0]=="@":
			if not has_gtf:
				# build hits dict
				if re.search("^\@SQ", szl):
					d = sam_header_parse(szl)
					dhits[d['SN']] = 0
					annot[d['SN']] = int(d['LN'])
			
			else:
				continue
		
		#
		# parse reads
		#aln = samaln_init(szl)
		#rname = aln['qname']
		aln = szl.strip().split("\t")
		# convert flag to int to avoid calls to int later in the code
		aln[SAM_FLAG] = int(aln[SAM_FLAG])
		
		rname = aln[SAM_QNAME]
		num_parsed += 1
		
		if rname != rlast:
			num_reads += 1
			if len(aln_buffer) > 0:
				num_aligned += 1
				# deal with alignment...
				rres = process_read(args, dhits, annot, aln_buffer)
				if rres[0] > 0:
					# update fragment mean and stdev
					frag_len_n += 1
					frag_len_mean = frag_len_mean0 + float(rres[1]-frag_len_mean0)/frag_len_n
					frag_len_stdev = frag_len_stdev0 + float(rres[1]-frag_len_mean0)*(rres[1]-frag_len_mean)
					frag_len_mean0 = frag_len_mean
					frag_len_stdev0 = frag_len_stdev
					
					num_assigned +=1
				
				if rres[0]==1:
					num_assigned_unique += 1
				elif rres[0]==2:
					num_assigned_multi += 1
			
			aln_buffer = []
		
		if (aln[SAM_FLAG] & SAM_UNAL) == 0:
			# don't bother with unaligned reads
			aln_buffer.append(aln)
			
		rlast = rname
		
		if (num_parsed % 1e6) == 0:
			message("Parsed {} fragments".format(num_reads))
	
	# close that fifo
	fin.close()
	
	if len(aln_buffer) > 0:
		num_aligned += 1
		# deal with final read
		rres = process_read(args, dhits, annot, aln_buffer)

		if rres[0] > 0:
			# update fragment mean and stdev
			frag_len_n += 1
			frag_len_mean = frag_len_mean0 + float(rres[1]-frag_len_mean0)/frag_len_n
			frag_len_stdev = frag_len_stdev0 + float(rres[1]-frag_len_mean0)*(rres[1]-frag_len_mean)
			frag_len_mean0 = frag_len_mean
			frag_len_stdev0 = frag_len_stdev
			
			num_assigned +=1
		
		if rres[0]==1:
			num_assigned_unique += 1
		elif rres[0]==2:
			num_assigned_multi += 1

	frag_len_stdev = sqrt(frag_len_stdev0/(frag_len_n-1))

	# increment for that last read that would be in the buffer (if aligned) when the 
	# main loop exits
	num_reads += 1	
	message("Parsed {} fragments".format(num_reads))

	sys.stderr.write("\n")
	sys.stderr.write("Mean fragment length: {:0.1f} +/- {:0.1f}\n\n".format(frag_len_mean, frag_len_stdev))
	sys.stderr.write("Parsed {} fragments of which:\n".format(num_reads))
	sys.stderr.write("    {} ({:0.1f}%) were aligned\n".format(num_aligned, num_aligned*100.0/num_reads))
	sys.stderr.write("\n")
	sys.stderr.write("Of {} fragments that passed filters:\n".format(num_aligned))
	sys.stderr.write("    {} ({:0.1f}%) were assigned to features with:\n".format(num_assigned, num_assigned*100.0/num_aligned))
	sys.stderr.write("      {} ({:0.1f}%) uniquely assigned to genes and\n".format(num_assigned_unique, num_assigned_unique*100.0/num_assigned))
	sys.stderr.write("      {} ({:0.1f}%) ambiguously assigned to genes\n".format(num_assigned_multi, num_assigned_multi*100.0/num_assigned))
	sys.stderr.write("\n")	

	# wait for fork
	#os.waitpid(childs[0], 0)

	try:
		os.unlink("{}.sam".format(ftmp))
	except:
		# haha
		message("haha")

	flist = ["{}.1.sai".format(ftmp), "{}.2.sai".format(ftmp)]
	for f in flist:
		if isfile(f):
			os.unlink(f)

	#
	# produce output
	#	
	
	if has_gtf:
		print "transcript_id\tgene_id\tgene_name\tchrom\tstrand\tlength\teff_length\thits\teff_hits\tfpkm\ttpm"
		dheader = { 'length':5, 'eff_length':6, 'eff_hits':8, 'fpkm':9 }
	else:
		print "transcript_id\tlength\teff_length\thits\teff_hits\tfpkm\ttpm"
		dheader = { 'length':1, 'eff_length':2, 'eff_hits':4, 'fpkm':5 }
	
	total_hits = 0
	total_fpkm = 0
	outlines = []
	for tid in sorted(dhits.keys()):
		#print "{}\t{}".format(tid, dhits[tid])
		lout = [tid]
		hits = dhits[tid]
		eff_hits = hits
		if has_gtf:
			length = annot[tid]['length']
			eff_length = annot[tid]['length']-frag_len_mean
			if eff_length <= 4:
				eff_length = 0
				hits = 0
				eff_hits = 0
			else:
				eff_hits = float(hits)*annot[tid]['length']/eff_length
			
			lout += [annot[tid]['gene_id'], annot[tid]['gene_name'], annot[tid]['rname'], annot[tid]['strand']]
		
		else:
			length = annot[tid]
			eff_length = annot[tid]-frag_len_mean 
			if eff_length <= 4:
				eff_length = 0
				hits = 0
				eff_hits = 0
			else:
				eff_hits = float(hits)*annot[tid]/eff_length
		
		lout += [length, eff_length, hits, eff_hits]
		total_hits += eff_hits
		outlines.append(lout)
		
		#print "\t".join(map(str, lout))
	
	# make fpkm
	nlines = len(outlines)
	for i in range(nlines):
		f = 0
		ll = outlines[i]
		if ll[dheader['eff_hits']] > 0 and ll[dheader['eff_length']] > 0:
			f = float(ll[dheader['eff_hits']])*1e9/(ll[dheader['eff_length']]*total_hits)
		
		total_fpkm += f
		outlines[i].append(f)
	
	# make tpm and print!
	for i in range(nlines):
		ll = outlines[i]
		tpm = 0
		if ll[-1] > 0:
			tpm = ll[dheader['fpkm']]*1e6/total_fpkm
		
		ll.append(tpm)
		for i in range(dheader['eff_length'], len(ll)):
			# format the floats to 4 decimals
			tmp = "{:0.4f}".format(ll[i])
			ll[i] = tmp
		print "\t".join(map(str, ll))
		

	# done!

	return 0

# the XA attribute....
# XA:Z:ENST00000429017,+113,100M,0;ENST00000548663,+191,100M,0;ENST00000450510,+1088,100M,0;ENST00000544999,+787,100M,0;

def parse_xa(aln):
	
	s0 = []
	saln = []
	xa0 = ""
	
	# find the XA tag if it exists
	for i in range(SAM_FIRST_ATTR, len(aln)):
		if re.search("^XA\:", aln[i]):
			# GOT IT
			xa0 = re.sub("^XA\:Z\:", "", aln[i])
			break
	
	xa0 = xa0.strip()

	if len(xa0) > 0:
		# string has length, continue
		s0 = xa0.split(";")
		# loop through split of the string
		for tmp in s0:
			if len(tmp) < 1:
				continue
				
			tmpHat = tmp.split(",")
			strand = tmpHat[XA_POS][0]
			tmpHat.append(strand)
			# strip the strand indicator off of the position
			tmp0 = re.sub("^[+-]", "", tmpHat[XA_POS])
			tmpHat[XA_POS] = int(tmp0)
			saln.append(tmpHat)
	
	return saln

def parse_xa_old(aln):
	
	s0 = []
	saln = []
	
	if "XA" in aln['attr']:
		# tag is present, get the string
		xa0 = aln['attr']['XA']
		if len(xa0) > 0:
			# string has length, continue
			s0 = xa0.split(";")
			# loop through split of the string
			for tmp in s0:
				if len(tmp) < 1:
					continue
					
				tmpHat = tmp.split(",")
				# second value has position and most importantly, strand
				atmp = { 
					'rname':tmpHat[0],
					'strand':tmpHat[1][0],
					 'pos':tmpHat[1],
					 'cigar':tmpHat[2]
				}
				atmp['pos'] = re.sub("^[+-]", "", atmp['pos'])
				saln.append(atmp)
	
	
	return saln
			
#
# this function handles a bundle of alignments all belonging to the same read 
def process_read(args, dhits, annot, buffer):
	
	# the flow:
	# pe or not?
	# if pe then 
	#   explode XA on each mate and append to mate
	#   pair mates
	#   count hits
	# else
	#   explode XA and assign hits
	
	pe = False
	tid_hit = set()
	gid_hit = set()
	flen = -1
	
	if (buffer[0][SAM_FLAG] & SAM_PAIRED) != 0:

		aleft = []
		aright = []
		# first of all there really should never be any more than 2 hits in here
		#if len(buffer) > 2:
		#	warning_message("More than two instances of a read name in PE data...how?")
			
		for i in range(len(buffer)):
			if (buffer[i][SAM_FLAG] & SAM_FIRST_MATE) != 0:
				aleft.append(buffer[i])
			else:
				aright.append(buffer[i])
		
		# left and right are sorted. first "pair" are the ones that were written
		# out as the primary alignment. we will also look into the secondary 
		# alignments. if primary alignment is a valid hit then we just have to 
		# check a couple things
		if len(aleft) == len(aright) and len(aleft)==1:
			# ok this is the expected result because BWA is supposed to only 
			# print a single primary alignment per pair with secondary hits
			# listed in the XA tag
						
			if aleft[0][SAM_RNAME] == aright[0][SAM_RNAME]:
				# good
				if args.fr_stranded:
					if sam_strand(aleft[0])=="-" and sam_strand(aright[0])=="+":
						# good!
						tid_hit.add(aleft[0][SAM_RNAME])
						# keep fragment length
						flen = abs(float(aleft[0][SAM_TLEN]))
						
				elif args.rf_stranded:
					if sam_strand(aleft[0])=="+" and sam_strand(aright[0])=="-":
						# good!
						tid_hit.add(aleft[0][SAM_RNAME])
						# keep fragment length
						flen = abs(float(aleft[0][SAM_TLEN]))
				else:
					# unstranded library. add the hit
					tid_hit.add(aleft[0][SAM_RNAME])
					# keep fragment length
					flen = abs(float(aleft[0][SAM_TLEN]))
						
			# 
			# expand the XA tags for each
			aleft_xa = parse_xa(aleft[0])
			aright_xa = parse_xa(aright[0])
			
			if len(aleft_xa) > 0 and len(aright_xa) > 0:
			
				#
				# now we have to pair these up and when we find a pair we can 
				# count the hit
				for i in range(len(aleft_xa)):
					aln_a = aleft_xa[i]
						
					for j in range(len(aright_xa)):
						aln_b = aright_xa[j]
						if aln_b is None:
							continue
						
						if aln_a[XA_RNAME]==aln_b[XA_RNAME]:
							# good!
							if args.fr_stranded:
								if (aln_a[XA_STRAND] != aln_b[XA_STRAND]) and (aln_a[XA_STRAND] == "-"):
									# good!
									tid_hit.add(aln_a[XA_RNAME])
									aleft_xa[i] = None
									aright_xa[j] = None
									continue
							elif args.rf_stranded:
								if (aln_a[XA_STRAND] != aln_b[XA_STRAND]) and (aln_a[XA_STRAND] == "+"):
									# good!
									tid_hit.add(aln_a[XA_RNAME])
									aleft_xa[i] = None
									aright_xa[j] = None
									continue
							else:
								if (aln_a[XA_STRAND] != aln_b[XA_STRAND]):
									# good!
									tid_hit.add(aln_a[XA_RNAME])
									aleft_xa[i] = None
									aright_xa[j] = None
									continue
							
			# done with that!
		elif len(aleft) > 1 or len(aright) > 1:
			# why?
			warning_message("More than two instances of a read name in PE data...how?")
		# else we probably are missing one of the mates
		
	else:
		# single-end data!
		
		for aln in buffer:
			aln_xa = parse_xa(aln)
			
			if args.fr_stranded:
				if sam_strand(aln) == "-":
					tid_hit.add(aln[SAM_RNAME])
					flen = sam_aligned_length(aln)
			elif args.rf_stranded:
				if sam_strand(aln) == "+":
					tid_hit.add(aln[SAM_RNAME])
					flen = sam_aligned_length(aln)
			else:
				# unstranded
				tid_hit.add(aln[SAM_RNAME])
				flen = sam_aligned_length(aln)
			
			for xaln in aln_xa:
				if args.fr_stranded:
					if xaln[XA_STRAND] == "-":
						tid_hit.add(xaln[XA_RNAME])
				elif args.rf_stranded:
					if xaln[XA_STRAND] == "+":
						tid_hit.add(xaln[XA_RNAME])
				else:
					# unstranded
					tid_hit.add(xaln[XA_RNAME])
			
	# 
	# now we have the 'tid_hit' set. change it to a list and get the corresponding gene id set
	tid_hit = list(tid_hit)
	if len(tid_hit) < 1:
		# no hits
		return [0, flen]
	
	#
	# look up the gene ids corresponding to these transcripts
	if args.G:
		for tid in tid_hit:
			if tid in annot:
				gid_hit.add(annot[tid]['gene_id'])
		
		gid_hit = list(gid_hit)
	
	elif args.L:
		for tid in tid_hit:
			if tid in annot:
				gid_hit.add(annot[tid]['attrs']['locus'])
		
		gid_hit = list(gid_hit)
		
	else:
		gid_hit = ['hello']
	
	
	
	#
	# now we have the ability to give this alignment a weight
	wi = 1
	ngenes = len(gid_hit)
	ntid = len(tid_hit)
	
	wi = 1.0/ntid * 1.0/ngenes**2 
	
	for tid in tid_hit:
		dhits[tid] += wi
	
	rres = 1
	if ngenes > 1:
		rres = 2
	
	return [rres, flen]

#
# launch subprocess to convert sai to sam. return the subprocess
# object so we can read from it's stdout in main
def bwa_do_sampe(args, stub, lsai):
	cmd = BWA + " sampe "
	cmd += "-a {} ".format(args.a)
	cmd += "-n {} ".format(args.N)
	cmd += "-f {}.sam ".format(stub)
	cmd += args.index
	cmd += " {} {} {} {} 2>sampe.log".format(lsai[0], lsai[1], args.reads[0], args.reads[1])
	#sys.stderr.write("CMD: {}\n".format(cmd))
	
	#p1 = sp.Popen(cmd.split(), stdout=sp.PIPE)
	#return p1
	return cmd

#
# launch subprocess to convert sai to sam. return the subprocess
# object so we can read from it's stdout in main
def bwa_do_samse(args, stub, lsai):
	cmd = BWA + " samse "
	cmd += "-n {} ".format(args.N)
	cmd += "-f {}.sam ".format(stub)
	cmd += args.index
	cmd += " {} {} 2>samse.log".format(lsai[0], args.reads[0])
	#sys.stderr.write("CMD: {}\n".format(cmd))
	
	#p1 = sp.Popen(cmd.split(), stdout=sp.PIPE)
	#return p1
	return cmd

#
# runs 'bwa aln' on input reads. writes to a .sai file
# named 'stub'.'mate'.sai
def bwa_do_alignment(read_queue, args, stub, lock):
	name = mp.current_process().name
	
	while True:
		task = read_queue.get()
		if task is None:
			read_queue.task_done()
			break
		
		reads = task[0]
		mate = task[1]
		
		message_mp("Aligning {}".format(reads), name, lock)
	
		cmd = BWA + " aln "
		cmd += "-n {} ".format(args.n)
		cmd += "-o {} ".format(args.o)
		cmd += "-i {} ".format(args.i)
		cmd += "-e {} ".format(args.e)
		cmd += "-t {} ".format(args.t)
		cmd += "-R {} ".format(args.R)
		cmd += "-l {} ".format(args.l)
		cmd += "-k {} ".format(args.k)
		cmd += "-f {}.{}.sai ".format(stub, mate)
		if args.I:
			cmd += "-I "
		
		cmd += args.index + " " + reads + " 2>bwa.log"
		message_mp("{}".format(cmd), name, lock)
		
		os.system(cmd)
		
		read_queue.task_done()
	
	return
	
	

def error_message(sz):
	sys.stderr.write("Error: {}\n".format(sz))

def warning_message(sz):
	sys.stderr.write("Warning: {}\n".format(sz))

def message(sz):
	sys.stderr.write("[bwa-quant] {}\n".format(sz))

def message_mp(sz, name, lock):
	lock.acquire()
	sys.stderr.write("[{}] {}\n".format(name, sz))
	lock.release()
	return

#==============================================================================
# SAM defs
#==============================================================================

# expecting that passed object is the list that was split from a single
# sam alignment line
def sam_aligned_length(laln):
	cigar = laln[SAM_CIGAR]
	r = re.findall("([0-9]+)[MD]", cigar)
	return sum(map(float, r))

def sam_strand(laln):
	flag = int(laln[SAM_FLAG])
	s = "+"
	if flag & SAM_REVERSED:
		s = "-"
	return s

def sam_are_mates(laln1, laln2):
	
	if laln1[SAM_RNEXT]==laln2[SAM_RNEXT]:
		if sam_strand(laln1) != sam_strand(laln2):
			if (int(laln1[SAM_PNEXT])==int(laln2[SAM_POS])) and (int(laln1[SAM_POS])==int(laln2[SAM_PNEXT])):
				return True
	
	return False 
	

def sam_header_parse(sz):
	ll = szl.strip().split("\t")
	d = {}
	d['type'] = re.sub("^@", "", ll[0])
	# parse additional fields
	for i in range(1, len(ll)):
		tmp = ll[i].split(":")
		d[tmp[0]] = tmp[1]
	
	return d
		

#
# parse a SAM record into a dict. if the read is aligned then the boundaries
# of the alignment are calculated as well as the aligned length
def samaln_init(sz):
	aln = {}

	tmp = sz.strip().split("\t")

	aln['qname'] = tmp[0]
	aln['flag'] = int(tmp[1])
	aln['rname'] = tmp[2]
	aln['pos'] = int(tmp[3])
	aln['mapq'] = int(tmp[4])
	aln['cigar'] = tmp[5]
	aln['rnext'] = tmp[6]
	aln['pnext'] = int(tmp[7])
	aln['tlen'] = int(tmp[8])
	aln['seq'] = tmp[9]
	aln['qual'] = tmp[10]

	aln['read_len'] = len(tmp[9])

	# find boundaries of the alignment as well as the aligned length which 
	# would be different from the read length in the event of soft-clipping
	if not samaln_unaligned(aln):
		samaln_calc_bounds(aln)

	aln['attr'] = {}

	# parse out attributes
	if len(tmp) > 11:
		for i in range(11, len(tmp)):
			tparts = tmp[i].split(":")
			value = tparts[2]
			if tparts[1] == "i":
				value = int(value)
			elif tparts[1] == "f":
				value = float(value)

			aln['attr'][tparts[0]] = value

	return aln

#
# return True if aln is the aligned mate of aln0
def samaln_is_mate(aln0, aln):

	if aln0['rnext'] == aln['rnext'] or aln0['rname'] == aln['rname']:
		# usually these are set to '=' when both alignments are in the same 
		# reference
		if aln0['pnext'] == aln['pos'] and aln0['pos'] == aln['pnext']:
			# positions match up
			return True

	return False




#
# functions to check the status of the alignment flag
def samaln_paired(aln):
	return True if (aln['flag'] & 0x1) != 0 else False

def samaln_properly_paired(aln):
	return True if (aln['flag'] & 0x1) != 0 else False

def samaln_unaligned(aln):
	return True if (aln['flag'] & 0x4) != 0 else False

def samaln_mate_unaligned(aln):
	return True if (aln['flag'] & 0x8) != 0 else False

def samaln_is_reversed(aln):
	return True if (aln['flag'] & 0x10) != 0 else False

def samaln_mate_reversed(aln):
	return True if (aln['flag'] & 0x20) != 0 else False

def samaln_is_first_mate(aln):
	return True if (aln['flag'] & 0x40) != 0 else False

def samaln_is_second_mate(aln):
	return True if (aln['flag'] & 0x80) != 0 else False

def samaln_is_secondary(aln):
	return True if (aln['flag'] & 0x100) != 0 else False

#
# compute the boundaries of an alignment
def samaln_calc_bounds(aln):
	aln['bounds'] = []

	left = aln['pos'] 
	right = aln['pos']

	op_len = re.split("[A-Z]", aln['cigar'])
	op_typ = re.split("[0-9]+", aln['cigar'])
	
	# drop the extra bit that comes along with the split
	op_len = map(int, op_len[0:(len(op_len)-1)])
	op_typ = op_typ[1:len(op_typ)]
	nop = len(op_typ)

	# loop through
	aln['aln_len'] = 0
	for i in range(nop):
		if op_typ[i] == "M" or op_typ[i] == "D":
			right += op_len[i]
			if op_typ[i] == "M":
				aln['aln_len'] += op_len[i]

		if op_typ[i] == "N":
			# start a new feature and pass the current on
			aln['bounds'].append([left, right-1])
			left = right + op_len[i]
			right = left

	# append feature
	aln['bounds'].append([left, right-1])

	return 0


#==============================================================================
# GTF defs
#==============================================================================

def gtf_parseline(sz):

	tmp = sz.strip().split("\t")

	grow = {
		'rname':tmp[0],
		'db':tmp[1],
		'type':tmp[2],
		'start':int(tmp[3]),
		'end':int(tmp[4]),
		'strand':tmp[6],
		'attrs':{}}

	# parse attributes
	fsplit = tmp[8].split("\"")
	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		grow['attrs'][key.strip()] = fsplit[i+1].strip()
		i += 2

	return grow


def gtf_transcript_id(grow):
	if "transcript_id" in grow['attrs']:
		return grow['attrs']['transcript_id']
	return None

def gtf_gene_id(grow):
	if "gene_id" in grow['attrs']:
		return grow['attrs']['gene_id']
	return None

def gtf_gene_name(grow):
	if "gene_name" in grow['attrs']:
		return grow['attrs']['gene_name']
	return None

def gtf_load(fname):
	#
	# load a gtf file into a few annotation tables

	dannot = {}
	dgid2tid = defaultdict(set)
	dgname2tid = defaultdict(set)

	fin = open(fname, "r")
	for szl in fin:
		grow = gtf_parseline(szl)

		tid = gtf_transcript_id(grow)
		if tid not in dannot:
			dannot[tid] = {
				'rname':grow['rname'],
				'strand':grow['strand'],
				'db':grow['db'],
				'gene_id':gtf_gene_id(grow),
				'gene_name':gtf_gene_name(grow),
				'num_exons':0,
				'start':grow['start'],
				'end':grow['end'],
				'length':0,
				'exons':[]}

		dannot[tid]['length'] += grow['end']-grow['start']+1
		dannot[tid]['exons'].append([grow['start'], grow['end']])
		dannot[tid]['num_exons'] += 1
		
		if grow['start'] < dannot[tid]['start']:
			dannot[tid]['start'] = grow['start']
		if grow['end'] > dannot[tid]['end']:
			dannot[tid]['end'] = grow['end']

		if dannot[tid]['gene_id'] is not None:
			dgid2tid[dannot[tid]['gene_id']].add(tid)

		if dannot[tid]['gene_name'] is not None:
			dgname2tid[dannot[tid]['gene_name']].add(tid)

	fin.close()

	#
	# pass through the annotation and sort the exons by position
	for tid in dannot.keys():
		dannot[tid]['exons'].sort(key=lambda x: x[0])

	return dannot, dgid2tid, dgname2tid

def gtf_num_exons(annot, tid):
	if tid not in annot:
		return -1
	
	return annot[tid]['num_exons']

def gtf_get_intron_chain(annot, tid):
	if tid not in annot:
		return None
	
	ints = []
	n = gtf_num_exons(annot, tid)
	if n < 2:
		return None

	elast = annot[tid]['exons'][0]
	for i in range(1, n):
		e = annot[tid]['exons'][i]
		ints.append([elast[1]+1, e[0]-1])
		elast = e
	
	return ints



def region_init(rname, start, end, tag=None, index=None):
	r = { 'rname':rname, 'start':start, 'end':end }
	
	if tag is not None:
		r['tag'] = tag

	if index is not None:
		r['index'] = index

	return r


#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(
	description="""Naive method of expression quantification using BWA aln to
                   align the reads to the transcriptome""")

parser.add_argument("index", type=str, 
	help="BWA aligner index")
parser.add_argument("reads", type=str, nargs="+", metavar="reads", 
	help="FASTQ or FASTA reads to align and quantify. If additional files then PE is assumed.")

# quantification options
parser.add_argument("-g", dest="g", type=str, action="store", default="", 
	help="GTF annotation. Automatically annotates output")
parser.add_argument("-L", dest="L", action="store_const", const=True, default=False,
	help="Enables multi-locus mapping ambiguity downweighting (with -g)")
parser.add_argument("-G", dest="G", action="store_const", const=True, default=False,
	help="Enables multi-geneid mapping ambiguity downweighting (with -g)")

parser.add_argument("--fr-stranded", action="store_const", const=True, default=False, 
	help="Library is stranded where first mate is from the reverse strand")
parser.add_argument("--rf-stranded", action="store_const", const=True, default=False, 
	help="Library is stranded where first mate is from the forward strand")

# bwa options
parser.add_argument("-n", dest="n", action="store", type=float, default=0.04,
	help="max #diff (int) or missing prob under 0.02 err rate (float) [0.04]")
parser.add_argument("-o", dest="o", action="store", type=int, default=1, 
	help="maximum number or fraction of gap opens [1]")
parser.add_argument("-e", action="store", type=int, default=-1, 
	help="maximum number or gap extensions, -1 for disabling long gaps [-1]")
parser.add_argument("-i", dest="i", action="store", type=int, default=5, 
	help="do not put an indel within INT bp towards the ends [5]")
parser.add_argument("-l", dest="l", action="store", type=int, default=24, 
	help="seed length [24]")
parser.add_argument("-k", dest="k", action="store", type=int, default=2, 
	help="maximum differences in the seed [2]")
parser.add_argument("-t", dest="t", action="store", type=int, default=1, 
	help="number of threads [1]")
parser.add_argument("-R", dest="R", action="store", type=int, default=30, 
	help="stop searching when there are >INT equally best hits [30]")
parser.add_argument("-I", dest="I", action="store_const", const=True, default=False, 
	help="the input is in the Illumina 1.3+ FASTQ-like format (base64)")

# bwa samse/sampe options:
parser.add_argument("-a", dest="a", action="store", type=int, default=500, 
	help="maximum insert size for PE [500]")
parser.add_argument("-N", dest="N", action="store", type=int, default=200, 
	help="maximum hits to output per alignment [200]")








args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

