# imports to avoid errors in eclipse
import sys, re, copy
import subprocess as sp


# class for *.tracking file rows
class TrackingRow(object):
	def __init__(self):
		self.tid = ""
		self.xloc = ""
		self.ref_gene = ""
		self.ref_tid = ""
		self.class_code = ""
		self.samples = []
		self.num_samples = 0
		self.present_samples = 0
		return None

	# 
	# this function parses the tracking file row into this class
	def parse(self, sz):
		sz = szl.strip()
		aln = sz.split("\t")
		n = len(aln)

		ref_split = aln[2].split("|")

		self.tid = aln[0]
		self.xloc = aln[1]
		self.ref_gene = ref_split[0]
		self.ref_tid = ref_split[1]
		self.class_code = aln[3]

		# copy the samples into a list
		self.samples = list(aln[4:n])
		# get number of samples
		self.num_samples = len(self.samples)
		# get count of samples that contain the feature
		for i in range(self.num_samples):
			if self.samples[i] != "-":
				self.present_samples += 1

		return 0


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


# class for junction objects
class Junction(object):
	def __init__(self, rname="", start=0, end=0, id="", index=0, 
				strand="+", count=0, tid=None, gname=None, gid=None):
		
		self.rname = rname
		self.start = int(start)
		self.end = int(end)
		self.id = id
		self.index = index
		self.strand = strand
		self.count = float(count)
		
		if tid is None:
			self.tid = []
		else:
			self.tid = list(tid)
		
		if gname is None:
			self.gname = []
		else:
			self.gname = list(gname)	

		if gid is None:
			self.gid = []
		else:
			self.gid = list(gid)
	
	# --
	# compare this to another junction object
	def compare(self, jobj):
		res = False
		
		if self.rname==jobj.rname and self.start==jobj.start and self.end==jobj.end:
			res = True
		
		return res
	
	#-- 
	# merge a second Junction into this one
	def merge(self, jobj):
		# not sure what to do with index when we merge these
		self.index = -1
		
		# update tid, gid and gname
		tmp = set(self.tid)
		tmp.update(jobj.tid)
		self.tid = list(tmp)
		
		tmp = set(self.gid)
		tmp.update(jobj.gid)
		self.gid = list(tmp)
		
		tmp = set(self.gname)
		tmp.update(jobj.gname)
		self.gname = list(tmp)
		
		return 0
	
	def set_id_from_position(self):
		self.id = "{}:{}-{}".format(self.rname, self.start, self.end)
		return 0
	
	# return string of the donor id
	def get_donor_id(self):
		did = "{}:{}:5p".format(self.rname, self.start)
		return did
	
	# return string of the acceptor id
	def get_acc_id(self):
		aid = "{}:{}:3p".format(self.rname, self.end)
		return aid
			

#--
# parse_gtf_to_junctions
# creates a list of Junction objects
# calling code should verify that the file exists and maybe throw this
# whole thing in a 'try'
def parse_gtf_to_junctions(f):
	
	ljuncs = []
	ltmp = []
	tid = ""
	last_tid = ""
	lrow = []
	aln = []
	szl = ""
	jidx = 0
	junc = None
	
	# open file
	fin = open(f, "r")
	
	# loop through gtf
	for szl in fin:
		aln = szl.strip().split("\t")
		attr = parse_gtf_attr(aln[8])
		tid = attr['transcript_id']
		
		if tid == last_tid:
			# make junction
			
			# keep index of junction
			jidx += 1
			# make a new Junction
			junc = Junction(rname=aln[0], start=int(lrow[4])+1, end=int(aln[3])-1, index=jidx, strand=aln[6])
			# make an id so that i don't have to do it later
			junc.set_id_from_position()
			# append in the id and the attributes
			junc.tid.append(attr['transcript_id'])
			junc.gid.append(attr['gene_id'])
			if "gene_name" in attr:
				junc.gname.append(attr['gene_name'])
				
			# add junction to the list - clearly there will be redundant junctions as we go
			ljuncs.append(copy.deepcopy(junc))
			
		else:
			jidx = 0
			
		lrow = list(aln)
		last_tid = tid 
	
	fin.close()
	
	return(ljuncs)

#
# junc_list_to_dict
# From a list of junctions to a dict keyed by the junction ids. matching junctions are 
# merged. annotation will propagate into junctions in this way.
def junc_list_to_dict(ljuncs):
	djuncs = {}
	n = len(ljuncs)
	
	# loop through
	for i in range(n):
		if ljuncs[i].id not in djuncs:
			djuncs[ljuncs[i].id] = copy.deepcopy(ljuncs[i])
		
		else:
			djuncs[ljuncs[i].id].merge(ljuncs[i])
	
	return djuncs

#
# parse a GTF into a dict of GtfRow objects
def parse_gtf(fname):
	# variables
	gtfdb = {}
	grow = None

	# open file and parse it
	fin = open(fname, "r")
	for szl in fin:

		grow = GtfRow()
		grow.parse(szl.strip())

		if grow.type != "exon":
			continue

		if grow.tid not in gtfdb:
			gtfdb[grow.tid] = []

		gtfdb[grow.tid].append(grow)

	fin.close()

	return gtfdb

def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	#
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

# --
# runcmd
# run a system level command in subprocess. optionally you can return the process.
# if the process isn't returned then the function waits for the process to finish
def runcmdp(cmd, returnProcess=False):
	sys.stderr.write("CMD: {}\n".format(cmd))
	p1 = sp.Popen(cmd.split())

	if returnProcess==True:
		return(p1)

	p1.wait()
	return(0)

#
# check for fai index for the supplied fasta
def fai_exists(fa):
	return(isfile("{}.fai".format(fa)))

#
# if the supplied fasta file doesn't have an fai index then this function
# calls samtools to build it
def build_fai(fa):
	cmd = "samtools faidx {}".format(fa)
	runcmd(cmd)
	return 0


#
# run a system level command. used for running alignment and samtools commands
def runcmd(cmd, verbose=True):
	if verbose:
		sys.stderr.write("CMD: {}\n".format(cmd))
	rres = os.system(cmd)
	return rres

#
# print message to stderr
def message(sz):
	sys.stderr.write("[<program name>] " + sz + "\n")


#--
# hash_region
# r is a list with [ref, start, end] for the region
# and bin is a binning integer for hashing.
# the region may fall into multiple bins. 
# returns a list of hashes
def hash_region(r, rbin):

	rbin = int(rbin)
	hstart = int(r[1])/rbin
	hend = int(r[2])/rbin

	hh = []
	for i in range(hstart, hend+1):
		hh.append("{}:{}".format(r[0], i))

	return(hh)

#--
# region_overlap
# compares regions a and b to see if they overlap. returns 
# a list: [0/1, overlap length]
def region_overlap(a, b):
	astart = float(a[1])
	aend = float(a[2])
	bstart = float(b[1])
	bend = float(b[2])
	ovl = 0
	olen = 0

	# check overlap
	if aend >= bstart and astart <= bend:
		ovl = 1
		olen1 = aend-bstart
		olen2 = bend-astart
		len1 = aend-astart
		len2 = bend-bstart
		# overlap length is the minimum of all of these lengths
		olen = min([olen1, olen2, len1, len2])

	return([ovl, olen])

#--
# build_point_lookup_hash
# lpoints is a list of, at minimum, two element lists that shall be hashed
def build_point_lookup_hash(lpoints, bin):

	pid = ""
	bin = int(bin)
	didx = {}
	
	for i in range(len(lpoints)):
		hpos = int(lpoints[i][1])/bin

		pid = "{}:{}".format(lpoints[i][0], hpos)

		if pid not in didx:
			didx[pid] = []

		# insert into the hash
		didx[pid].append(list(lpoints[i]))


	return(didx)

#--
# build_range_lookup_hash
# lpoints is a list of, at minimum, two element lists that shall be hashed
def build_range_lookup_hash(lranges, bin):

	rid = []
	bin = int(bin)
	didx = {}

	for i in range(len(lranges)):
		hh = hash_region(lranges[i], bin)

		for j in range(len(hh)):
			rid = hh[j]

			if rid not in didx:
				didx[rid] = []

			# insert into the hash
			didx[rid].append(list(lranges[i]))

	return(didx)


#--
# find_point_hash_hits
# look up a point or region in a point hash
def find_region_hash_hits(a, h, bin):

	# hash a, a region
	hh = hash_region(a, bin)
	# check for hits
	hids = {}

	for i in range(len(hh)):
		if hh[i] in h:
			# maybe, get list of elements at this node
			lcand = h[hh[i]]
			for j in range(len(lcand)):
				# overlap?
				rres = region_overlap(a, lcand[j])
				if rres[0]==1:
					# keep this hit
					hid = range_to_id(lcand[j])
					hids[hid] = 0

	return(hids.keys())


#--
# find_point_hash_hits
# look up hits to points for a range, a
def find_point_hash_hits(a, h, bin):

	# hash a, a region
	hh = hash_region(a, bin)
	# check for hits
	hids = {}

	for i in range(len(hh)):
		if hh[i] in h:
			# maybe, get list of elements at this node
			lcand = h[hh[i]]
			for j in range(len(lcand)):
				# overlap?
				rres = region_overlap(a, lcand[j])
				if rres[0]==1:
					# keep this hit
					hid = point_to_id(lcand[j])
					hids[hid] = 0

	return(hids.keys())



def range_to_id(r):
	rid = "{}:{}-{}".format(r[0], r[1], r[2])
	return(rid)

def point_to_id(r):
	pid = "{}:{}".format(r[0], r[1])
	return(pid)

#--
# parse_bam_alignment
# parse a bam alignment into a dict
def parse_bam_alignment(sz):
	# vars
	aln = []
	fields = ["qname", "flag", "rname", "pos", "mapq", "cigar", 
		"rnext", "pnext", "tlen", "seq", "qual", "ext"]
	ext_idx = 11
	daln = {}

	# explode
	aln = sz.strip().split("\t")
	nfields = len(aln)

	for i in range(0, ext_idx):
		daln[fields[i]] = aln[i]

	daln["ext"] = list(aln[ext_idx:nfields])

	return(daln)

#--
# parse_sam_header_row
# parses a sam header line into a list
def parse_sam_header_row(sz):
	aln = sz.strip().split("\t")
	dh = {}
	x = aln[0]
	dh["type"] = x[1:len(x)]
	ll = []
	for i in range(1,len(aln)):
		tmp = aln[i].split(":")
		ll.append(tmp[1])

	dh["data"] = ll
	return(dh)


#--
# explode_cigar
# explode cigar string into operations and lengths
# note: uses re
def explode_cigar(cigar):
	
	# using findall we can blow this thing up all in one shot
	res = []
	if cigar != "*":
		res = re.findall("([0-9]+)([MIDNSHPX=]+)", cigar)
	
	return(res)

def aln_is_spliced(daln):
	rres = False
	rr = re.search("[0-9]+N", daln["cigar"])
	if rr:
		rres = True
	return rres

#--
# aln_get_left_right
# get the left and right most coordinates (adding soft clips back in)
# for the alignment. expect a dict input.
def aln_get_left_right(daln):

	res = [-1, -1]
	left0 = 0
	left = 0
	right = 0

	if daln['cigar'] != "*" and aln_is_aln(daln):
		left0 = int(daln['pos'])
		left = left0
		right = left0

		# adjust left for soft-clipping?
		rres = re.search("^([0-9]+)S", daln["cigar"])
		if rres:
			left -= int(rres.group(1))

		# get ops that change position
		rres = re.findall("([0-9]+)([MDN])", daln["cigar"])
		# find right
		for i in range(len(rres)):
			right += int(rres[i][0])

		# adjust right for soft-clipping?
		rres = re.search("([0-9]+)S$", daln["cigar"])
		if rres:
			right += int(rres.group(1))

		res = [left, right]

	return res


#--
# aln_get_mapped_regions
# get the left and right most coordinates (adding soft clips back in)
# for the alignment. expect a dict input.
def aln_get_mapped_regions(daln):

	mrout = []
	left0 = 0
	left = 0
	right = 0

	if daln['cigar'] != "*" and aln_is_aln(daln):
		left = int(daln['pos'])

		# get ops that change position
		rres = re.findall("([0-9]+)([MDN=])", daln["cigar"])
		# find "M" ops or "=" ops
		for i in range(len(rres)):
			if rres[i][1] == "M" or rres[i][1] == "=":
				# make a region
				right = left+int(rres[i][0])-1
				mrout.append([daln["rname"], left, right])

			# move left position
			left += int(rres[i][0])

	return mrout

#--
# get a dict of True/False for all of the possible flags in an 
# alignment.
def aln_check_flags(daln):
	rres = {}
	flag = int(daln["flag"])

	rres["pe"] = (flag & 0x1)
	rres["proper"] = (flag & 0x2)
	rres["unmapped"] = (flag & 0x4)
	rres["mate_unmapped"] = (flag & 0x8)
	rres["reversed"] = (flag & 0x10)
	rres["mate_reversed"] = (flag & 0x20)
	rres["first_mate"] = (flag & 0x40)
	rres["second_mate"] = (flag & 0x80)
	rres["secondary"] = (flag & 0x100)
	rres["filtered"] = (flag & 0x200)
	rres["duplicate"] = (flag & 0x400)
	rres["supaln"] = (flag & 0x800)

	return(rres)


#--
# aln_is_aln
# return true if the alignment is aligned
def aln_is_aln(daln):
	res = True
	if int(daln["flag"]) & 0x4:
		res = False
	return(res)
