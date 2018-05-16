#!/usr/bin/python
#==============================================================================
# mt-file-process.py
#
# Shawn Driscoll
# 20170505
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Example of processing an input file with multiple queue workers. This example
# parses a GTF file and bundles transcript ids by gene name 
#==============================================================================

import sys, argparse, math, re
from os.path import isfile, expanduser
from os import system
import multiprocessing
import itertools

# use this as a way to combine all dicts
from collections import defaultdict

# use for timing blocks of code
import timeit
import time

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

def mains(args):
	
	# 
	t0 = timeit.default_timer()
	
	grow = None
	results = defaultdict(set)
	num_lines = 0
	
	with open(args.gtf, "r") as fin:
		iters = itertools.chain(fin)
		for l in iters:
			grow = GtfRow()
			grow.parse(l)
			results[grow.gid].update([grow.tid])
			num_lines += 1
	
	t1 = timeit.default_timer()
	
	part1 = t1-t0
	sys.stderr.write("parsed {} lines in {:0.4f} seconds\n".format(num_lines, part1))

	gid_list = results.keys()
	ngid = len(gid_list)
	ntid = 0
	for gid in gid_list:
		ntid += len(set(results[gid]))	
		
	sys.stderr.write("number of genes: {}; number of transcripts: {}\n".format(ngid, ntid))

	return 0	

def main(args):

	t0 = timeit.default_timer()

	# create a pool
	pool = multiprocessing.Pool(args.threads)
	
	num_lines = 0
	max_lines = 2**13
	buffer = []
	results = []
	dres = defaultdict(list)
	
	# open input file
	fin = open(args.gtf, "r")
	for szl in fin:
		
		num_lines += 1		
		buffer.append(szl)
		
		if len(buffer) >= max_lines:
			results = pool.map(generic_pool_worker, buffer)
			join_dict(dres, results)
			buffer = []

	fin.close()

	if len(buffer) > 0:
		results = pool.map(generic_pool_worker, buffer)
		join_dict(dres, results)
		buffer = []
	
	t1 = timeit.default_timer()
	
	pool.close()
	pool.join()
	
	part1 = t1-t0
	sys.stderr.write("parsed {} lines in {:0.4f} seconds\n".format(num_lines, part1))

	gid_list = dres.keys()
	ngid = len(gid_list)
	ntid = 0
	for gid in gid_list:
		ntid += len(set(dres[gid]))	
		
	sys.stderr.write("number of genes: {}; number of transcripts: {}\n".format(ngid, ntid))
	
	# done!
	return 0
		

# receives a list of dicts	
def join_dict(dout, ld):
	
	for d in ld:
		for k in d.keys():
			dout[k] += d[k]
	
	return 0
		

# receives a single line from file and deals with it
def generic_pool_worker(item):
	grow = GtfRow()
	grow.parse(item)
	
	return { grow.gid: [ grow.tid ] }


def main_foobar(args):
		

	t0 = timeit.default_timer()
	
	tasks = multiprocessing.JoinableQueue()
	results = multiprocessing.Queue()
	lk = multiprocessing.Lock()
	pool = []
	num_lines = 0
	dres = defaultdict(set)
	chunk_count = 0
	chunk_lines = 2**8
	
	last_line = None
	
	#
	# create worker pool
	for i in range(args.threads):
		p = multiprocessing.Process(target=generic_worker, args=(tasks, results))
		p.daemon = True
		p.start()
		pool.append(p)
	
	fin = open(args.gtf, "r")
	for szl in fin:
		num_lines += 1
		chunk_count += 1
			
#		if (num_lines % 1000) == 0:
#			sys.stderr.write("parsed {} lines\n".format(num_lines))
		tasks.put(szl)
		last_line = szl
		
		if chunk_count >= chunk_lines:
			# process results queue
#			for p in pool:
#				tasks.put(None)

#			for p in pool:
#				p.join()
			
			sys.stderr.write("max lines hit, processing results queue...\n")
			sys.stderr.write("joining tasks\n")
			
			# this only blocks until all tasks removed by 'get' have been 
			# finished with 'task_done()' calls.
			tasks.join()
			
			results.put(None)

			sys.stderr.write("looping over results...\n")
			while True:			
				lout = results.get()
				
				if lout is None:
					break
					
				dres[lout[0]].add(lout[1])
				
			#time.sleep(10)
			chunk_count = 0
		
#			for p in pool:
#				p.start()

			# continue...
	
	fin.close()

	# pass None to the queue to signal the end of the data	
	for i in range(args.threads):
		sys.stderr.write("sending poison pill...\n")
		tasks.put(None)

	sys.stderr.write("parsed {} lines\n".format(num_lines))
	
#	sys.stderr.write("joining tasks queue\n")
#	tasks.join()

	sys.stderr.write("joining processes\n")
	for p in pool:
#		print p
		p.join()

	sys.stderr.write("finished, wrapping up\n")

	t1 = timeit.default_timer()
	
	part1 = t1-t0
	sys.stderr.write("parsed {} lines in {:0.4f} seconds\n".format(num_lines, part1))
	
	# should be done, right?
	results.put(None)
	
	while True:
		lout = results.get()
		if lout is None:
			break
		
		dres[lout[0]].add(lout[1])
	
	gid_list = dres.keys()
	ngid = len(gid_list)
	ntid = 0
	for gid in gid_list:
		ntid += len(set(dres[gid]))	
		
	sys.stderr.write("number of genes: {}; number of transcripts: {}\n".format(ngid, ntid))
	
	return 0
			

def generic_worker(task_queue, result_queue):
	
	name = multiprocessing.current_process().name
	
	while True:
		item = task_queue.get()
		
		if item is None:
			task_queue.task_done()
			break
		
#		sys.stderr.write("a\n")
		tsk = Task(item)
#		sys.stderr.write("b\n")
		
		rres = tsk()
#		sys.stderr.write("c\n")
		
		result_queue.put(rres)
#		sys.stderr.write("d\n")
		
		task_queue.task_done()
#		sys.stderr.write("e\n")
	
	sys.stderr.write("exiting process {}\n".format(name))
	
	return
		

def main_manager(args):
	
	t0 = timeit.default_timer()
	
	#
	# "manual" map pool approach. 
	manager = multiprocessing.Manager()
	results = manager.list()
	work = manager.Queue(args.threads)
	num_lines = 0
	
	#
	# start workers
	pool = []
	for i in range(args.threads):
		p = multiprocessing.Process(target=pool_worker, args=(work, results))
		p.start()
		pool.append(p)
	
	with open(args.gtf, "r") as fin:
		# reads the file, 'fin', and then inserts 'None' for as many 
		# processes we have running
		iters = itertools.chain(fin, (None,)*args.threads)
		
		# loops through the file...
		for l in iters:
			work.put(l)
			num_lines += 1
		
	for p in pool:
		p.join()

	t1 = timeit.default_timer()
	
	# done!
	dout = defaultdict(set)
	for r in results:
		dout[r[0]].update([r[1]])

	t2 = timeit.default_timer()
	
	part1 = t1-t0
	part2 = t2-t1
	
	sys.stderr.write("parsed {} lines in {:0.4f} seconds\npart2 {:0.4f}\n".format(num_lines, part1, part2))

	gid_list = dout.keys()
	ngid = len(gid_list)
	ntid = 0
	for gid in gid_list:
		ntid += len(set(dout[gid]))	
		
	sys.stderr.write("number of genes: {}; number of transcripts: {}\n".format(ngid, ntid))

	
def pool_worker(in_queue, dout):
	while True:
		item = in_queue.get()
		
		if item is None:
			return
		
#		time.sleep(0.01)
		grow = GtfRow()
		grow.parse(item)
		
		# update list
		dout.append([grow.gid, grow.tid])

def main2(args):

	# variables
	tasks = multiprocessing.JoinableQueue()
	results = multiprocessing.JoinableQueue()
	processes = []

	#
	# create processes and fire them up
	processes = [DaFunk(tasks, results) for i in range(args.threads)]
	for p in processes:
		p.start()
	
	rproc = DaResults(results)
	rproc.start()

	fin = open(args.gtf, "r")
	lset = []
	aln = []
	gid = ""
	gidLast = ""
	dout = defaultdict(set)
	lk = multiprocessing.Lock()
	
	tid = ""
	dtid = defaultdict(int)
	
	max_buff = 2**16 - args.threads*2

	t0 = timeit.default_timer()
	idx = 0
	for szl in fin:
		# push lines straight into the queue. the nature of this example
		# gives us the opportunity to process fewer than all of the lines 
		# from the file since once a transcript id is encountered we can 
		# skip all other instances of it

		tasks.put(Task(szl))
				
		# parse tid
#		r = re.search("transcript_id \"([^\"]+)\"", szl)
#		if r:
#			tid = r.group(1)
#			if tid not in dtid:
#				dtid[tid] += 1
#				# process this line
#				hq.q.put(szl)

	fin.close()

	# kill processes
	for i in range(args.threads):
		tasks.put(None)

	# join
	tasks.join()

	loop_time = timeit.default_timer() - t0
	sys.stderr.write("finished parsing file in {:0.4}\n".format(loop_time))

	results.put(None)
	
	results.join()
	
	print rproc.result
	
#	while True:
#		l = results.get()
#		if l is None:
#			# finished!
#			break
#		dout[l[0]].update([l[1]])
#	
#	print len(dout.keys())
#	nid = 0
#	for k in dout.keys():
#		nid += len(dout[k])
#	print nid
	
	return 0


#==============================================================================
# class
#==============================================================================

#
# This class is based on a process and adds in the in and out queues as well 
# as a lock just in case. the 'run' function is triggered when you call 
# the '.start()' function on instances of this object
class DaFunk(multiprocessing.Process):
	
	def __init__(self, task_queue, result_queue, lock):
		multiprocessing.Process.__init__(self)
		self.task_queue = task_queue
		self.result_queue = result_queue
		self.lock = lock
	
	def run(self):
		proc_name = self.name
		while True:
			item = self.task_queue.get()
			
			# break out if we hit a None
			if item is None:
				self.task_queue.task_done()
				break
			
			# process the item
			rres = item()
			
			self.task_queue.task_done()
			#self.lock.acquire()
			self.result_queue.put(rres)
			#self.lock.release()
		
		return


class DaResults(multiprocessing.Process):
	
	def __init__(self, result_queue):
		multiprocessing.Process.__init__(self)
		self.result_queue = result_queue
		self.result = defaultdict(set)
	
	def run(self):
		proc_name = self.name
		while True:
			item = self.result_queue.get()
			
			# break out if we hit a None
			if item is None:
				self.result_queue.task_done()
				break
			
			# process the item
			self.result[item[0]].update([item[1]])
			print self.result
			
			self.result_queue.task_done()
		
		return

#
# task class can be defined to do all of the work necessary for a single task.
# the '__call__' function makes it so you can run the object like a function
# with no parameters
class Task(object):
	def __init__(self, item):
		self.item = item
	
	def __call__(self):
		grow = GtfRow()
		grow.parse(self.item)
		return [grow.gid, grow.tid]

#
# file processing queue. distributes file processing tasks to multiple workers
# and builds output in a dict. one output per worker. this class should also
# provide a function to merge the multiple outputs into a single output
class FPQueue(object):

	# init
	def __init__(self, threads):
		# class members
		self.threads = threads
		self.q = Queue()
		self.qdone = Queue()
		self.l = Lock()
		self.processes = []
		# empty dict for the thread outputs
		self.results = defaultdict(set)

		#
		# start threads
		for i in range(threads):
			# start processes
			p = Process(target=self.worker)
			p.daemon = True
			p.start()
			# keep them handy in case we have trouble
			self.processes.append(p)

		# good

	def process_done(self):
		for k in iter(self.qdone.get, None):
			self.results[k[0]].update([k[1]])
		
		return 0
				

	# stop processes by pushing in a 'None' object
	def complete_job(self):
		# join queue
#		self.q.join()		
		# send stop commands

		for p in self.processes:
			self.q.put(None)
		
		for p in self.processes:
			p.join()

		# join again
#		self.q.join()
		
		for i in range(self.threads):
			self.qdone.put(None)
			
		self.process_done()
		
		return self.results
				
	# main worker
	def worker(self):
		# get name of this thread
		name = current_process().name
		# make sure this thread has a dict entry
#		if name not in self.results:
#			self.results[name] = {}
#		if name not in dres:
#			dres[name] = {}

		for item in iter(self.q.get, None):
			self.process_item(item, name)
#			self.q.task_done()
		
#		self.q.task_done()
		
		return True

	def process_item(self, item, tname):
		# deal with the item. in this case the item is a single line from
		# a gtf file that we need to parse and exract the gene id and transcript id

		grow = GtfRow()
		grow.parse(item)

#		self.l.acquire()
#		if grow.gid not in self.results[tname]:
#			self.results[tname][grow.gid] = []

		# update gene id with the current transcript id
#		self.results[tname][grow.gid].append(grow.tid)
#		self.l.release()
		
#		self.l.acquire()
#		self.results[grow.gid].update([grow.tid])
		#print grow.gid, self.results[grow.gid]
#		self.l.release()

		self.qdone.put([grow.gid, grow.tid])
		
		return 0

	def merge_results(self):

		# dout is the final merged dict
		tnames = self.results.keys()

		# The results, which are all dicts in this case, can be merged
		# either blindly or carefully. blindly if there are not going to be
		# common key ids in the separate dicts and carefully if there will be

		# merge dicts via response found here: 
		# http://stackoverflow.com/questions/1781571/how-to-concatenate-two-dictionaries-to-create-a-new-one-in-python
		if False:
			dout = dict(self.results[tnames[0]], **self.results[tnames[1]])
			for tname in tnames[2:len(tnames)]:
				dout.update(self.results[tname])

		if False:
			dout = {}
			for tname in self.results.keys():
				dout.update(self.results[tname])

		# normal loop
		if True:

			dout = defaultdict(set)

			# combine dicts into a list
			ld = []
			for t in self.results.keys():
				ld.append(self.results[t])

			for d in ld:
				for key, value in d.iteritems():
					dout[key].update(value)

		# done!
		return dout


#
# class to hold GTF row objects. 
class GtfRow(object):
	def __init__(self, rname="", db="", type="", start=0, end=0, strand="."):
		self.rname = rname
		self.db = db
		self.type = type
		self.start = start
		self.end = end
		self.strand = strand
		
		self.tid = ""
		self.gid = ""
		
		self.attrs = {}
	
	def __str__(self):
		lout = self.tolist()
		return "\t".join(map(str, lout))

	# parse info in GTF row, sz, into this object
	def parse(self, sz):
		aln = sz.strip().split("\t")
		self.rname = aln[0]
		self.db = aln[1]
		self.type = aln[2]
		self.start = int(aln[3])
		self.end = int(aln[4])
		self.strand = aln[6]
		
		# parse attributes from field 8
		fsplit = aln[8].split("\"")
		n = len(fsplit)-1
		i = 0
		while i < n:
			key = re.sub(';','',fsplit[i])
			self.attrs[key.strip()] = fsplit[i+1].strip()
			i += 2

		self.tid = self.attrs['transcript_id']
		self.gid = self.attrs['gene_id']
		
		return 0

	#--
	# tolist
	# return a list version of this with elements in place of the columns 
	# of a GTF
	def tolist(self):
		ltmp = [self.rname, self.db, self.type, 
			int(self.start), int(self.end), ".", self.strand, "."]
		
		szattr = "transcript_id \"{}\"; gene_id \"{}\";".format(self.tid, self.gid)
		akey = self.attrs.keys()
		if len(akey) > 0:
			# append additional attributes to the szattr string
			for aid in akey:
				szattr += " {} \"{}\";".format(aid, self.attrs[aid])
		
		ltmp.append(szattr)
		
		return(ltmp)


#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Example multithreaded file processing.")
parser.add_argument('gtf', type=str, help="Input GTF")
parser.add_argument('-t', '--threads', type=int, default=1, 
	help="Threads for parallel processing.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

