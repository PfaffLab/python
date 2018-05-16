#!/usr/bin/env python
#==============================================================================
# snpa.py
#
# shawn driscoll
# 20120731
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parse output of post-snp pipeline (filter-snps) to find the exact codon
# of the mutation. Reference codon and mutant codon are returned as well as 
# their amino translation
#==============================================================================

import sys,argparse,re
import HTSeq as hts
import math

#==============================================================================
# functions
#==============================================================================

##
# revComp
# Returns the reverse compliment of a DNA string, sz
def revComp(sz):
	dr = {"A":"T","C":"G","T":"A","G":"C"}
	n = len(sz)
	i = 0
	
	# convert sequence to a list, make a copy
	lsz = list(sz)
	tmp = ""
	
	# loop backward through the sequence, compliment and insert
	# into copy starting from first position
	while n >= i:
		n -= 1
		tmp = dr[lsz[n]]
		lsz[n] = dr[lsz[i]]
		lsz[i] = tmp
		i += 1

	return("".join(lsz))


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError:
		return False

	fin.close()
	return True

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
# make a fast lookup hash of the GTF
#
def gtf_lookup_hash(fname, ftype):

	lk = {}

	#- open file -#

	fin = open(fname, "r")

	#- parse file and build hash

	for szl in fin:
		ll = szl.strip().split("\t")
		
		if ll[2] == ftype:
			kid = hash_pos(ll[0], ll[3])

			if kid not in lk:
				lk[kid] = []

			lk[kid].append(list(ll))

	fin.close()

	return(lk)

def hash_pos(rname, pos):
	mod = 16000
	bin = int(pos)/mod
	kid = rname + str(bin)
	return kid

# 
# look up a position in the hash, return list of intersections
#
def lookup_position(lk, lpos):

	result = []

	query_id = hash_pos(lpos[0], lpos[1])

	if query_id in lk:
		# check for intersections with features in this bucket
		ll = lk[query_id]
		for i in range(len(ll)):
			lf = ll[i]
			if int(lf[3]) <= int(lpos[1]) and int(lf[4]) >= int(lpos[1]):
				# append entire GTF record so that the main function can 
				# have access to all feature information
				result.append(list(lf))

	return result


#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Look up codons from SNPs")
parser.add_argument('snps',type=str,help="SNP results file processed by filter-snps")
parser.add_argument('gtf',type=str,help="GTF gene reference annotation that includes CDS regions")
parser.add_argument('fa',type=str,help="FASTA genome reference file. File has to have an index created by fasta-index.py in the same directory")
parser.add_argument("-a",action="store_const",dest="b_printAll",const=True,default=False,help="Print all rows from input file (default: CDS rows only)")
parser.add_argument("--changed-only",action="store_const",dest="b_changedOnly",const=True,default=False,help="Print only rows where the amino acid was altered (default: False)")
parser.add_argument("--min-hits",dest="n_minHits",action="store",default=2,type=int,help="Minimum high quality hits for a SNP to be included in the analysis (and output) (default: 2)")

args = parser.parse_args()

# check for files
try:
	fin = open(args.snps,"r")
	fin.close()
except IOError as e:
	print "Error: unable to find/open input file: " + args.snps
	sys.exit()

try:
	fin = open(args.gtf,"r")
	fin.close()
except IOError as e:
	print "Error: unable to find/open input file: " + args.gtf
	sys.exit()

try:
	fin = open(args.fa,"r")
	fin.close()
except IOError as e:
	print "Error: unable to find/open input file: " + args.fa
	sys.exit()

try:
	fin = open(args.fa + ".fai","r")
	fin.close()
except IOError as e:
	print "Error: unable to find/open FASTA index file: " + args.fa + ".fai"
	sys.exit()


#==============================================================================
# variables
#==============================================================================

# genitic code reference table
sz_gCode = "/Users/pfafflab/projects/genome/gcode.txt"

dGtf = {}
dFai = {}
szl = ""
ll = []
szTid = ""
fin = None
sz_seq = ""
sz_tmp = ""
l_tmp = []
l_tindex = []
l_t = []
b_hasCds = False
b_print = True
n_ti = 0
n_codonPosition = 0
sz_codon = ""
dGcode = {}
l_genes = []
l_cdsi = []
n_hits = 0
n_hqhits = 0
l_snpBases = []

has_gene_names = False
tid_to_gene_name = {}

b_debug = False

#==============================================================================
# main script
#==============================================================================

##
# read in genitic code index
try:
	fin = open(sz_gCode,"r")
	for szl in fin:
		szl = szl.strip()
		ll = szl.split("\t")
		dGcode[ll[0]] = ll[1:]

	fin.close()
except IOError,e:
	sys.stderr.write("{0} ({1})\n",e.strerror,sz_gCode)
	sys.exit()

##
# read in fasta index
sys.stderr.write("> parsing reference index...\n")
fin = open(args.fa + ".fai","r")
for szl in fin:
	szl = szl.strip()
	ll = szl.split("\t")
	dFai[ll[0]] = ll[2:]

fin.close()

###
# parse GTF
sys.stderr.write("> parsing " + args.gtf + " for CDS records...\n")
fin = hts.GFF_Reader(args.gtf)
for feature in fin:
	if feature.type == "CDS":
		szTid = feature.attr['transcript_id']

		if szTid not in dGtf:
			dGtf[szTid] = []

		# can't count of GTF files to be sorted correctly so we have to be 
		# careful while inserting features into each transcript's list
		if len(dGtf[szTid]) > 0:
			i = 0
			while feature.iv.start > dGtf[szTid][i].start:
				i += 1
				if i >= len(dGtf[szTid]):
					break

			if i > len(dGtf[szTid]):
				dGtf[szTid].append(feature.iv)
			else:
				dGtf[szTid].insert(i, feature.iv)

		else:
			dGtf[szTid].append(feature.iv)
			
	elif feature.type == "exon":
		szTid = feature.attr['transcript_id']

		if "gene_name" in feature.attr:
			has_gene_names = True
			if szTid not in tid_to_gene_name:
				tid_to_gene_name[szTid] = set([])

			tid_to_gene_name[szTid].update([feature.attr['gene_name']])


sys.stderr.write("> parsing SNP data and locating codons...\n")

# open snps file
fin = open(args.snps,"r")
# open fasta file
ffa = open(args.fa,"r")

###
# loop through snps file and handle each one
for szl in fin:
	szl = szl.strip()
	ll = szl.split("\t")

	# figure out if this SNP hit a CDS, track each CDS position so that
	# we can also track which gene names correspond to CDS hits
	l_cdsi = []
	l_t = ll[8].split(",")
	b_print = True
	b_hasCds = False
	for i in range(len(l_t)):
		if l_t[i] == "CDS":
			b_hasCds = True
			#n_ti = i
			l_cdsi.append(i)

	# modify 8th column (index 7) and pick out the read count values
	n_hits = 0
	n_hqhits = None
	m = re.search("DP\=([0-9]+)",ll[7])
	if m:
		n_hits = int(m.group(1))

	m = re.search("DP4\=([0-9\,]+)",ll[7])
	if m:
		l_t = m.group(1).split(",")
		n_hqhits = sum(map(int,l_t[2:]))

	ll[7] = "\t".join(map(str,[n_hits,n_hqhits]))


	# bail out if this row has two (or more) different non-ref base values
	#if len(ll[4]) > 1:
	#	b_hasCds = False

	# compare total depth at this snp to minimum depth
	if n_hits < args.n_minHits:
		b_print = False
		b_hasCds = False

	if b_hasCds:
		##
		# This SNP hit a CDS...continue

		# modify the gene column to include only gene names that have a CDS for this
		# hit (in the case of loci with multiple gene names)
		l_t = ll[9].split(",")
		l_genes = []
		for i in range(len(l_cdsi)):
			l_genes.append(l_t[l_cdsi[i]])

		ll[9] = ",".join(list(set(l_genes)))

		# take first CDS index for the transcript id we'll use for analysis
		n_ti = l_cdsi[0]
		l_t = ll[10].split(",")
		szTid = l_t[n_ti]

		# initalize some stuff
		l_tindex = []  # used to build an index of the transcript 
		sz_seq = ""    # storage for the sequence of this transcript
		n_tpos = 1     # start position of transcript base at first feature
		n_snpPos = int(ll[1])  # snp position in the genome


		# fetch sequence for the transcript and build the transcript index, l_tindex
		# as we go
		read_length = 0
		for fiv in dGtf[szTid]:
			##
			# fetch position of the chromosome from the fasta index file as well as the 
			# line width and actual line width (including end of line junk)
			n_chr = int(dFai[fiv.chrom][0])
			n_l = int(dFai[fiv.chrom][1])
			n_la = int(dFai[fiv.chrom][2])
			# append index row with genomic coordinates for this feature as well as 
			# converted feature coordinates
			l_tindex.append([fiv.start,fiv.end,n_tpos,(n_tpos+(fiv.end-fiv.start)-1)])
			n_tpos += fiv.end-fiv.start

			# n_chr is the position to seek to in the FASTA file, now find the 
			# offset to the start of this feature
			n_estart = int(fiv.start)
			# integer division of the start point by the line length and multiplying it by 
			# the actual line length gives the position of the line within the fasta file
			# where this feature's sequence starts
			n_estartAdj = n_estart/n_l*n_la
			# there will be some remainder bases (offset into the line)
			n_estartRem = n_estart - (n_estart/n_l)*n_l
			# calculate feature length
			n_elen = int(fiv.end)-n_estart
			n_elenAdj = n_elen/n_l*n_la
			n_elenRem = n_elen - (n_elen/n_l)*n_l

			read_length += n_elen

			# seek file pointer
			ffa.seek(n_chr+n_estartAdj+n_estartRem)
			# read sequence with some extra to account for end of line junk even
			# though it seems like we did this already
			sz_tmp = ffa.read(n_elenAdj+n_elenRem+max(1,int(math.ceil(n_elen/float(n_l)))))
#			print sz_tmp
			# split result by newlines
			l_tmp = sz_tmp.split("\n")
			# merge back together without newlines
			sz_tmp = "".join(l_tmp)
#			print "[NOTE] read_length=%d" % len(sz_tmp)
			# slice out the full feature length and append to final sequence
			sz_seq += sz_tmp[0:n_elen]
#			print "[NOTE] n_elen=%d; cumulative=%d; actual=%d" % (n_elen, read_length, len(sz_seq))

#		if b_debug:
#			print sz_seq

		# convert sequence to all upper case
		sz_seq = sz_seq.upper()

		# check seq length compared to what l_tindex says 
		if len(sz_seq) != l_tindex[-1][3]:
			print l_tindex
			print szl
			sys.stderr.write("[snpa] Error: sequencing length does not match expected length:\n")
			sys.stderr.write("       %s,%d found %d\n" % (szTid, l_tindex[-1][3], len(sz_seq)))
			sys.exit()

		##
		# we need to convert the snp's genomic coordinate into an intra CDS
		# coordinate so first we find which CDS chunk the SNP is within
		n_ti = -1
		for i in range(len(l_tindex)):
			if n_snpPos > l_tindex[i][0] and n_snpPos <= l_tindex[i][1]:
				n_ti = i
				break

		if n_ti < 0:
			sys.stderr.write("Error: snp doesn't hit this feature\n")
			sys.stderr.write("Line: %s\n" % szl)
			sys.exit()

		# convert position by finding it's offset from the start of the containing feature
		# and then adding that to the intra-CDS start of the same feature (minus 1 for some reason)
		n_snpFeatureOffset = n_snpPos-l_tindex[n_ti][0]-1
		n_snpCdsOffset = l_tindex[n_ti][2]+n_snpFeatureOffset-1
		# determin the positions relative position with the codon. this value will be 0, 1 or 2
		# for start, middle and end of the codon
		n_codonPosition = n_snpCdsOffset % 3
		
		# list of mutant bases (usually only one)
		l_snpBases = ll[4].split("/")
		# start mutant codon list
		l_mutantCodon = ["" for i in range(len(l_snpBases))]

		#sz_mutant = ""

		if b_debug:
			#print sz_seq
			print len(sz_seq)
			print l_tindex
			print szTid
			print fiv.chrom,n_snpPos,n_snpCdsOffset,n_codonPosition

		# pull together the reference and mutant version of the codon for each SNP base in l_snpBases
		ni = 0
		for ni in range(len(l_snpBases)):
			sz_mbase = l_snpBases[ni]

			if n_codonPosition == 0:
				sz_codon = sz_seq[n_snpCdsOffset] + sz_seq[n_snpCdsOffset+1] + sz_seq[n_snpCdsOffset+2]
				l_mutantCodon[ni] = sz_mbase + sz_seq[n_snpCdsOffset+1] + sz_seq[n_snpCdsOffset+2]
			elif n_codonPosition == 1:
				sz_codon = sz_seq[n_snpCdsOffset-1] + sz_seq[n_snpCdsOffset] + sz_seq[n_snpCdsOffset+1]
				l_mutantCodon[ni] = sz_seq[n_snpCdsOffset-1] + sz_mbase + sz_seq[n_snpCdsOffset+1]
			elif n_codonPosition == 2:
				sz_codon = sz_seq[n_snpCdsOffset-2] + sz_seq[n_snpCdsOffset-1] + sz_seq[n_snpCdsOffset]
				l_mutantCodon[ni] = sz_seq[n_snpCdsOffset-2] + sz_seq[n_snpCdsOffset-1] + sz_mbase

			if b_debug:
				print sz_codon,l_mutantCodon

		# reverse compliment if the strand is negative
		if fiv.strand == "-":
			sz_codon = revComp(sz_codon)
			for ni in range(len(l_mutantCodon)):
				l_mutantCodon[ni] = revComp(l_mutantCodon[ni])

		#if b_debug:	
			#print sz_codon,sz_mutant
			#print ",".join(dGcode[sz_codon]), ",".join(dGcode[sz_mutant])

		#print sz_codon,sz_mutant,dGcode[sz_codon],dGcode[sz_mutant]
		#print "\t".join(ll) + "\t" + "\t".join([sz_codon,dGcode[sz_codon][1],sz_mutant,dGcode[sz_mutant][1]]) + "\t" + fiv.strand
		# build string for mutant genitic code output
		l_mcodes = []
		for ni in range(len(l_mutantCodon)):
			try:
				l_mcodes.append(dGcode[l_mutantCodon[ni]])
			except KeyError:
				sys.stderr.write("[snpa] Error: key value does not exist (%s)\n" % l_mutantCodon[ni])
				sys.stderr.write("       Line: %s\n" % szl)
				sys.exit()

		# generate string version of the genitic code data
		sz_mutant = ""
		for li in l_mcodes:
			if len(sz_mutant) > 0:
				sz_mutant += "|"

			sz_mutant += ",".join(li)

		l_wtcode = dGcode[sz_codon]
		
		# insert hit ratio

		ll.insert(8, str(float(n_hqhits)/float(n_hits)))
		# simplify columns for exon, gene_id and transcript_id
		ll[9] = "CDS"
		if has_gene_names:
			# replace gene_id column with gene name
			# split the transcript ids
			ltmp = ll[11].split(",")
			gtmp = set([])
			for i in range(len(ltmp)):
				gtmp.update(tid_to_gene_name[ltmp[i]])

			ll[10] = ",".join(list(gtmp))
		else:
			ltmp = list(set(ll[10].split(",")))
			ll[10] = ",".join(ltmp)
		ltmp = list(set(ll[11].split(",")))
		ll[11] = ",".join(ltmp)

		if args.b_changedOnly:
			b_changed = False
			for ni in range(len(l_mutantCodon)):
				if l_wtcode[0] != l_mcodes[ni][0]:
					b_changed = True
			
			if b_changed:
				sys.stdout.write("\t".join(ll))
				sys.stdout.write("\t" + sz_codon + "\t" + ",".join(l_wtcode))
				sys.stdout.write("\t" + "|".join(l_mutantCodon) + "\t" + sz_mutant)
				sys.stdout.write("\n")

			#if dGcode[sz_codon][0] != dGcode[sz_mutant][0]:
			#	sys.stdout.write("\t".join(ll))
			#	sys.stdout.write("\t" + sz_codon + "\t" + ",".join(dGcode[sz_codon]))
			#	sys.stdout.write("\t" + sz_mutant + "\t" + ",".join(dGcode[sz_mutant]))
			#	sys.stdout.write("\n")
		else:
			sys.stdout.write("\t".join(ll))
			#sys.stdout.write("\t" + sz_codon + "\t" + ",".join(dGcode[sz_codon]))
			#sys.stdout.write("\t" + sz_mutant + "\t" + ",".join(dGcode[sz_mutant]))
			sys.stdout.write("\t" + sz_codon + "\t" + ",".join(l_wtcode))
			sys.stdout.write("\t" + "|".join(l_mutantCodon) + "\t" + sz_mutant)
			sys.stdout.write("\n")
	else:
		ll.insert(8, str(float(n_hqhits)/float(n_hits)))
		# simplify columns for exon, gene_id and transcript_id
		ll[9] = "exon"
		ltmp = list(set(ll[10].split(",")))
		ll[10] = ",".join(ltmp)
		ltmp = list(set(ll[11].split(",")))
		ll[11] = ",".join(ltmp)

		if args.b_printAll and b_print:
			sys.stdout.write("\t".join(ll))
			sys.stdout.write("\t".join(["-","-"]))
			sys.stdout.write("\n")

ffa.close()
