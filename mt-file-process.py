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
from Queue import Queue
import threading

from multiprocessing import cpu_count

# use this as a way to combine all dicts
from collections import defaultdict

# use for timing blocks of code
import timeit

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
	hq = FPQueue(threads=args.threads)

	fin = open(args.gtf, "r")
	lset = []
	aln = []
	gid = ""
	gidLast = ""

	tid = ""
	dtid = defaultdict(int)

	t0 = timeit.default_timer()
	
	for szl in fin:

		# push lines straight into the queue. the nature of this example
		# gives us the opportunity to process fewer than all of the lines 
		# from the file since once a transcript id is encountered we can 
		# skip all other instances of it
		
		# parse tid
		r = re.search("transcript_id \"([^\"]+)\"", szl)
		if r:
			tid = r.group(1)
			if tid not in dtid:
				dtid[tid] += 1
				# process this line
				hq.q.put(szl)

	fin.close()

	loop_time = timeit.default_timer() - t0
	
	# wait for queue to finishe
	t0 = timeit.default_timer()
	hq.q.join()
	join_time = timeit.default_timer() - t0


	for k in hq.results.keys():
		sys.stderr.write(k)
		sys.stderr.write("\t")
		sys.stderr.write(str(len(hq.results[k].keys())))
		sys.stderr.write("\n")


	# pull results together
	t0 = timeit.default_timer()
	rres = hq.merge_results()
	merge_time = timeit.default_timer()-t0

	sys.stderr.write("\nloop time {:f}\n".format(loop_time))
	sys.stderr.write("join time {:f}\n".format(join_time))
	sys.stderr.write("merge time {:f}\n\n".format(merge_time))

	# write results to stdout
	for k in sorted(rres.keys()):
		print k + "\t" + ",".join(rres[k])	

	return 0


#==============================================================================
# class
#==============================================================================

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
		# empty dict for the thread outputs
		self.results = {}

		#
		# start threads
		for i in range(threads):
			# create thread running the class's worker function
			t = threading.Thread(target=self.worker)
			t.daemon = True
			t.start()

		# good

	# main worker
	def worker(self):
		# get name of this thread
		name = threading.currentThread().getName()
		# make sure this thread has a dict entry
		if name not in self.results:
			self.results[name] = {}

		# main loop
		while True:
			# fetch item
			item = self.q.get()
			# deal with item...
			self.process_item(item, name)
			# when we're finished...
			self.q.task_done()

	def process_item(self, item, tname):
		# deal with the item. in this case the item is a single line from
		# a gtf file that we need to parse and exract the gene id and transcript id

		grow = GtfRow()
		grow.parse(item)

		if grow.gid not in self.results[tname]:
			self.results[tname][grow.gid] = []

		# update gene id with the current transcript id
		self.results[tname][grow.gid].append(grow.tid)

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

