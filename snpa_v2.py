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
# import HTSeq as hts
import math
import os

# -- 
# GLOBALS
# --

_HASH_MOD = 16000

#==============================================================================
# functions
#==============================================================================

##
# revComp
# Returns the reverse compliment of a DNA string, sz
def revComp(sz):
	dr = {"A":"T","C":"G","T":"A","G":"C", "N":"N", "-":"-"}
	
	lsz = list(sz)
	n = len(lsz)

	if n > 0:
		
		# reverse the list
		lsz.reverse()
		
		# replace nc's
		for i in range(n):
			lsz[i] = dr[lsz[i]]

	
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
#
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
# parse_gtf
# parses gtf CDS into a hash keyed by 'transcript_id'
#
def parse_gtf(fname):
	
	gtf_cds = {}
	gtf_exons = {}
	szl = ""
	ll = []
	attr = {}
	i = 0
	n = 0
	
	n_cds = 0
	n_exons = 0
	
	fin = open(fname, "r")
	
	for szl in fin:
		ll = szl.strip().split("\t")
		
		if ll[2] == "CDS":
			attr = parse_gtf_attr(ll[8])
			tid = attr['transcript_id']
			
			if tid not in gtf_cds:
				gtf_cds[tid] = []
				n_cds += 1
			
			if len(gtf_cds[tid]) == 0:
				gtf_cds[tid].append([ll[0], int(ll[3])-1, int(ll[4])])
			else:
				# find insertion to the list so the exons are sorted properly
				ltmp = [ll[0], int(ll[3])-1, int(ll[4])]
				i = 0
				while ltmp[1] > gtf_cds[tid][i][1]:
					i += 1
					if i >= len(gtf_cds[tid]):
						break
				
				gtf_cds[tid].insert(i, list(ltmp))
		elif ll[2] == "exon":
			attr = parse_gtf_attr(ll[8])
			tid = attr['transcript_id']
			
			if tid not in gtf_exons:
				gtf_exons[tid] = []
				n_exons += 1
			
			if len(gtf_exons[tid]) == 0:
				gtf_exons[tid].append(list(ll))
			else:
				# find insertion to the list so the exons are sorted properly
				i = 0
				while int(ll[3]) > int(gtf_exons[tid][i][3]):
					i += 1
					if i >= len(gtf_exons[tid]):
						break
				
				gtf_exons[tid].insert(i, list(ll))
	
	fin.close()
	sys.stderr.write("[parse_gtf] parsed %d exon and %d CDS features\n" % (n_exons, n_cds))
	
	return gtf_cds, gtf_exons

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
			# kid = hash_pos(ll[0], ll[3])
			kid_list = hash_feature([ll[0], ll[3], ll[4]])

			for kid in kid_list:
				if kid not in lk:
					lk[kid] = []
	
				lk[kid].append(list(ll))

	fin.close()

	return(lk)

def hash_pos(rname, pos):
	bucket = int(pos)/_HASH_MOD
	kid = rname + str(bucket)
	return kid

def hash_feature(lf):
	# --
	# considers the possibility that a feature would hash into several bins,
	# which totally happens with introns because sometimes they are long
	# -- 

	mstart = int(lf[1])/_HASH_MOD
	mend = int(lf[2])/_HASH_MOD

	kid_list = []

	kid_list.append(lf[0]+str(mstart))
	
	i = mstart+1
	while i <= mend:
		kid_list.append(lf[0]+str(i))
		i += 1

	return kid_list
# 
# look up a position in the hash, return list of intersections - designed for SNPs
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

def parse_vcf_info(sz):
	info = {}
	
	ltmp = sz.split(";")
	for i in range(len(ltmp)):
		if len(ltmp[i]) > 0:
			ltmp2 = ltmp[i].split("=")
			info[ltmp2[0]] = ltmp2[1]
	
	return info

def list_sum(l):
	l_n = map(float, l)
	return(sum(l_n))

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Look up SNPs relative to a GTF annotation returning codons and mutant codons for SNPs within coding sequences of genes (requires CDS features in GTF). Default behavior is to print out only CDS hits.")
parser.add_argument('snps',type=str,help="VCF SNP call output (expected from GATK)")
parser.add_argument('gtf',type=str,help="GTF gene reference annotation that includes CDS regions")
parser.add_argument('fa',type=str,help="FASTA genome reference file. File has to have an index created by fasta-index.py in the same directory")
parser.add_argument("-a",action="store_const",dest="b_printAll",const=True,default=False,help="Print all rows from input file (default: CDS rows only)")
parser.add_argument("--genes-only", action="store_const", dest="b_genes_only", const=True, default=False, 
				help="Print out all hits in genes ignoring intronic hits (default: off)")
parser.add_argument("--changed-only",action="store_const",dest="b_changedOnly",const=True,default=False,
				help="Print only SNPs with codon changes in CDS of a gene (default: off)")
parser.add_argument("--min-hits",dest="n_minHits",action="store",default=2,type=int,
				help="Minimum high quality hits for a SNP to be included in the analysis (and output) (default: 2)")


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
sz_gCode = os.getenv("HOME") + "/projects/genome/gcode.txt"

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
l_vcfInfo = []
num_raw_cols = 0
acounts = []
num_samples = 0
sample_names = []
sample_genotypes = []
sample_depths = []
sample_altRatios = []
genotypes = {"0/0": "negative", "0/1": "het", "1/1": "mutant", "./.":"null"}



# look up hashes for CDS and for exon features from the GTF
lk_exon = {}
lk_cds = {}

has_gene_names = False
tid_to_gene_name = {}

#b_debug = True
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
	sys.stderr.write("{0} ({1})\n".format(e.strerror,sz_gCode))
	sys.exit()

##
# read in fasta index
sys.stderr.write("[main] parsing reference index...\n")
fin = open(args.fa + ".fai","r")
for szl in fin:
	szl = szl.strip()
	ll = szl.split("\t")
	dFai[ll[0]] = ll[2:]

fin.close()

###
# parse GTF
sys.stderr.write("[main] parsing %s\n" % args.gtf)
dGtf, eGtf = parse_gtf(args.gtf)

#if False:
#	fin = hts.GFF_Reader(args.gtf)
#	for feature in fin:
#		if feature.type == "CDS":
#			szTid = feature.attr['transcript_id']
#	
#			if szTid not in dGtf:
#				dGtf[szTid] = []
#	
#			# can't count of GTF files to be sorted correctly so we have to be 
#			# careful while inserting features into each transcript's list
#			if len(dGtf[szTid]) > 0:
#				i = 0
#				while feature.iv.start > dGtf[szTid][i].start:
#					i += 1
#					if i >= len(dGtf[szTid]):
#						break
#	
#				if i > len(dGtf[szTid]):
#					dGtf[szTid].append(feature.iv)
#				else:
#					dGtf[szTid].insert(i, feature.iv)
#	
#			else:
#				dGtf[szTid].append(feature.iv)
#				
#		elif feature.type == "exon":
#			szTid = feature.attr['transcript_id']
#	
#			if "gene_name" in feature.attr:
#				has_gene_names = True
#				if szTid not in tid_to_gene_name:
#					tid_to_gene_name[szTid] = set([])
#	
#				tid_to_gene_name[szTid].update([feature.attr['gene_name']])
#
#

#
# - make lookup hashes for exons and CDS features -
#

sys.stderr.write("[main] building fast lookup tables of exons and CDS features...\n")
lk_exon = gtf_lookup_hash(args.gtf, "exon")
lk_cds = gtf_lookup_hash(args.gtf, "CDS")

sys.stderr.write("[main] parsing SNP data and locating codons...\n")

# open snps file
fin = open(args.snps,"r")
# open fasta file
ffa = open(args.fa,"r")

###
# loop through snps file and handle each one
for szl in fin:

	if re.search("^\#CHROM", szl):
		# header row
		ll = szl.strip().split("\t")
		num_samples = len(ll)-9
		sys.stdout.write("\t".join(ll[0:8] + ["feature", "exon_type", "gene_name", "transcript_id", "ref_codon", "ref_amino", "mut_codon", "mut_amino"]))
		for i in range(num_samples):
			sys.stdout.write("\t" + ll[9+i] + "_gtype")
		for i in range(num_samples):
			sys.stdout.write("\t" + ll[9+i] + "_depth")
		for i in range(num_samples):
			sys.stdout.write("\t" + ll[9+i] + "_altDepth")

		sys.stdout.write("\n")


		continue

	elif szl[0] == "#":
		continue

	ll = szl.strip().split("\t")
	lout = list(ll[0:8] + ["-" for i in range(8)])
	
	if len(ll[3]) > 1:
		# this is not a SNP but maybe a multi-base substitution
		continue

	szl = szl.strip()
	ll = szl.split("\t")
	
	if num_raw_cols == 0:
		num_raw_cols = len(ll) 

	num_samples = len(ll)-9

	# break down the genotype and depth ratios per sample (assuming GATK output)

#	sample_genotypes = []
#	sample_depths = []
#	sample_altRatios = []
#	for i in range(num_samples):
#		tmp = ll[9+i]
#		tmp_split = tmp.split(":")
#		sample_genotypes.append(genotypes[tmp_split[0]])
#
#		sample_depths.append(0)
#		sample_altRatios.append(0)
#
#		if tmp_split[0] != "./.":
#			if tmp_split[1] != ".":
#				# get hit counts
#				tmp2 = tmp_split[1].split(",")
#				rcount = int(tmp2[0])
#				acount = int(tmp2[1])
#				if rcount > 0 or acount > 0:
#					sample_depths[-1] = int(tmp2[0])+int(tmp2[1])
#					sample_altRatios[-1] = float(tmp2[1])/(int(tmp2[0])+int(tmp2[1]))
#	
#	nhits = 0
#	for d in sample_depths:
#		nhits += d

	# trim the sample results down to just genotype and count depths
	gtypes = ["-" for i in range(num_samples)]
	depths = [0 for i in range(num_samples)]
	alt_depths = [0 for i in range(num_samples)]
	for i in range(num_samples):
		tmp = ll[9+i]
		tmp_split = tmp.split(":")

		gtypes[i] = tmp_split[0]
		if len(tmp_split) > 1:
			if len(tmp_split[1]) > 0:
				if tmp_split[1] != ".":
					dtmp = tmp_split[1].split(",")
					try:
						depths[i] = int(list_sum(dtmp))
						alt_depths[i] = int(list_sum(dtmp[1:len(dtmp)]))
					except ValueError as e:
						sys.stderr.write("ValueError when summing depths: ")
						sys.stderr.write(str(e))
						sys.stderr.write("; " + str(tmp) + "\n")



	lout += gtypes
	lout += map(str,depths)
	lout += map(str,alt_depths)
	

	# figure out if this SNP hit a CDS, track each CDS position so that
	# we can also track which gene names correspond to CDS hits
#	l_cdsi = []
#	l_t = ll[8].split(",")

	b_print = True
	b_hasCds = False
	b_inGene = False
	
	s_cds_genes = set([])
	s_cds_tids = set([])
	s_exon_genes = set([])
	s_exon_tids = set([])
	s_exon_types = set([])
	
	cds_hits = lookup_position(lk_cds, [ll[0], ll[1]])
	exon_hits = lookup_position(lk_exon, [ll[0], ll[1]])
	
	b_hasCds = len(cds_hits) > 0
	b_inGene = len(exon_hits) > 0
	
	# pull gene and transcript names for this snp
	if b_hasCds:
		for i in range(len(cds_hits)):
			attr = parse_gtf_attr(cds_hits[i][8])
			if 'gene_name' in attr:
				s_cds_genes.update([attr['gene_name']])
			else:
				s_cds_genes.update([attr['gene_id']])
			
			s_cds_tids.update([attr['transcript_id']])
	
	if b_inGene:
		for i in range(len(exon_hits)):
			attr = parse_gtf_attr(exon_hits[i][8])
			if 'gene_name' in attr:
				s_exon_genes.update([attr['gene_name']])
			else:
				s_exon_genes.update([attr['gene_id']])
			
			s_exon_tids.update([attr['transcript_id']])
			
			# find this exon in the eGtf table
			tid = attr['transcript_id']
			tl = eGtf[tid]
			for j in range(len(tl)):
				if exon_hits[i][3] == tl[j][3] and exon_hits[i][4] == tl[j][4]:
					# got it
					if j == 0 and tl[j][6] == "+":
						if len(tl) == 1:
							s_exon_types.update(["se"])
						else:
							s_exon_types.update(["5p"])
					elif j == 0 and tl[j][6] == "-":
						if len(tl) == 1:
							s_exon_types.update(["se"])
						else:
							s_exon_types.update(["3p"])
					elif j == len(tl)-1 and tl[j][6] == "+":
						s_exon_types.update(["3p"])
					elif j == len(tl)-1 and tl[j][6] == "-":
						s_exon_types.update(["5p"])
					else:
						s_exon_types.update(["internal"])
					
					break
					
						
	
#	if len(cds_hits) > 0:
#		b_hasCds = True
#		# get feature names 
#		for i in range(len(lhits)):
#			attr = parse_gtf_attr(lhits[i][8])	
	
#	for i in range(len(l_t)):
#		if l_t[i] == "CDS":
#			b_hasCds = True
#			#n_ti = i
#			l_cdsi.append(i)

	# parse vcf info
#	l_vcfInfo = parse_vcf_info(ll[7])
#	n_hits = 1
#	if "DP" in l_vcfInfo:
#		n_hits = int(l_vcfInfo["DP"])
#	else:
#		sys.stderr.write("[main] warning: no 'DP' value @ snp: %s\n" % szl)

	# see near start of this loop where we figure out read depths
	
	
		
	# append gene names and transcript names
	if b_hasCds or b_inGene:
		lout[9] = ",".join(list(s_exon_types))
		lout[10] = ",".join(list(s_cds_genes))
		lout[11] = ",".join(list(s_cds_tids))
		
	if b_hasCds:
		lout[8] = "CDS"
	elif b_inGene:
		lout[8] = "exon"


	# establish whether or not we need to print this out
	b_print = True

	# have a look at the observed alternate alleals. check for INDELs which have alt strings 
	# longer than a single base

	a_bases = ll[4].split(",")
	a_bases_mod = []
	if len(a_bases) > 1:
		for i in range(len(a_bases)):
			if len(a_bases[i]) == 1:
				a_bases_mod.append(a_bases[i])
	elif len(a_bases[0]) == 1:
		a_bases_mod.append(a_bases[0])

	if len(a_bases_mod) > 0:
		ll[4] = ",".join(a_bases_mod)
	else:
		b_hasCds = False
		b_print = False

		
	if args.b_changedOnly and not b_hasCds:
		b_print = False 
	elif not b_inGene and args.b_genes_only:
		# if not in a gene and we only want gene hits then this is a no-printer
		b_print= False
	

	if b_hasCds:
		##
		# This SNP hit a CDS...continue

		# modify the gene column to include only gene names that have a CDS for this
		# hit (in the case of loci with multiple gene names)
#		l_t = ll[9].split(",")
#		l_genes = []
#		for i in range(len(l_cdsi)):
#			l_genes.append(l_t[l_cdsi[i]])
#
#		ll[9] = ",".join(list(set(l_genes)))

#		l_genes = []
#		for i in range(len(cds_hits)):
#			attr = parse_gtf_attr(cds_hits[i][8])
#			# snag the transcript id so we can fetch this snp's codon
#			if i == 0:
#				szTid = attr['transcript_id']
#			if "gene_name" in attr:
#				l_genes.append(attr['gene_name'])
#			else:
#				l_genes.append(attr['gene_id'])

#		# take first CDS index for the transcript id we'll use for analysis
		attr = parse_gtf_attr(cds_hits[0][8])
		szTid = attr['transcript_id']
		t_strand = cds_hits[0][6]
#		n_ti = l_cdsi[0]
#		l_t = ll[10].split(",")
#		szTid = l_t[n_ti]

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
			# n_chr = int(dFai[fiv.chrom][0])
			n_chr = int(dFai[fiv[0]][0])
			# n_l = int(dFai[fiv.chrom][1])
			n_l = int(dFai[fiv[0]][1])
			# n_la = int(dFai[fiv.chrom][2])
			n_la = int(dFai[fiv[0]][2])
			# append index row with genomic coordinates for this feature as well as 
			# converted feature coordinates
			# l_tindex.append([fiv.start,fiv.end,n_tpos,(n_tpos+(fiv.end-fiv.start)-1)])
			l_tindex.append([fiv[1],fiv[2],n_tpos,(n_tpos+(fiv[2]-fiv[1])-1)])
			n_tpos += (fiv[2]-fiv[1])

			# n_chr is the position to seek to in the FASTA file, now find the 
			# offset to the start of this feature
			# n_estart = int(fiv.start)
			n_estart = fiv[1]
			# integer division of the start point by the line length and multiplying it by 
			# the actual line length gives the position of the line within the fasta file
			# where this feature's sequence starts
			n_estartAdj = n_estart/n_l*n_la
			# there will be some remainder bases (offset into the line)
			n_estartRem = n_estart - (n_estart/n_l)*n_l
			# calculate feature length
			# n_elen = int(fiv.end)-n_estart
			n_elen = fiv[2]-fiv[1]
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
		l_snpBases = ll[4].split(",")
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
			
			try:
				if n_codonPosition == 0:
					sz_codon = sz_seq[n_snpCdsOffset] + sz_seq[n_snpCdsOffset+1] + sz_seq[n_snpCdsOffset+2]
					l_mutantCodon[ni] = sz_mbase + sz_seq[n_snpCdsOffset+1] + sz_seq[n_snpCdsOffset+2]
				elif n_codonPosition == 1:
					sz_codon = sz_seq[n_snpCdsOffset-1] + sz_seq[n_snpCdsOffset] + sz_seq[n_snpCdsOffset+1]
					l_mutantCodon[ni] = sz_seq[n_snpCdsOffset-1] + sz_mbase + sz_seq[n_snpCdsOffset+1]
				elif n_codonPosition == 2:
					sz_codon = sz_seq[n_snpCdsOffset-2] + sz_seq[n_snpCdsOffset-1] + sz_seq[n_snpCdsOffset]
					l_mutantCodon[ni] = sz_seq[n_snpCdsOffset-2] + sz_seq[n_snpCdsOffset-1] + sz_mbase
			except IndexError as e:
				#print len(sz_seq), n_snpCdsOffset
				#print e
				sz_codon = "---"


			if b_debug:
				print sz_codon,l_mutantCodon

		# reverse compliment if the strand is negative
		# if fiv.strand == "-":
		if t_strand == "-":
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
#				try:
			if l_mutantCodon[ni] in dGcode:
				l_mcodes.append(dGcode[l_mutantCodon[ni]])
			else:
				l_mcodes.append(["unk"])
#				except KeyError:
#					sys.stderr.write("[main] Error: key value does not exist (%s)\n" % l_mutantCodon[ni])
#					sys.stderr.write("       Line: %s\n" % szl)
#					sys.exit()

		# generate string version of the genitic code data
		sz_mutant = ""
		for li in l_mcodes:
			if len(sz_mutant) > 0:
				sz_mutant += "|"

			sz_mutant += ",".join(li)

		try:
			l_wtcode = dGcode[sz_codon]
		except KeyError as e:
			#print e
			l_wtcode = "-"
		
		# insert hit ratio

#		ll.insert(8, str(float(n_hqhits)/float(n_hits)))
#		# simplify columns for exon, gene_id and transcript_id
#		ll[9] = "CDS"
#		if has_gene_names:
#			# replace gene_id column with gene name
#			# split the transcript ids
#			ltmp = ll[11].split(",")
#			gtmp = set([])
#			for i in range(len(ltmp)):
#				gtmp.update(tid_to_gene_name[ltmp[i]])
#
#			ll[10] = ",".join(list(gtmp))
#		else:
#			ltmp = list(set(ll[10].split(",")))
#			ll[10] = ",".join(ltmp)
#		ltmp = list(set(ll[11].split(",")))
#		ll[11] = ",".join(ltmp)

		if args.b_changedOnly:
			b_changed = False
			for ni in range(len(l_mutantCodon)):
				if len(l_wtcode)==3 and len(l_mcodes[ni])==3:
					if l_wtcode[2] != l_mcodes[ni][2]:
						b_changed = True
				else:
					if l_wtcode[0] != l_mcodes[ni][0]:
						b_changed = True

			
			if not b_changed:
				b_print = False

		lout[12] = sz_codon
		lout[13] = ",".join(l_wtcode)
		lout[14] = ",".join(l_mutantCodon)
		lout[15] = sz_mutant
		
#		ll += [sz_codon, ",".join(l_wtcode), "|".join(l_mutantCodon), sz_mutant]
#			
#			if b_changed:
#				sys.stdout.write("\t".join(ll))
#				sys.stdout.write("\t" + sz_codon + "\t" + ",".join(l_wtcode))
#				sys.stdout.write("\t" + "|".join(l_mutantCodon) + "\t" + sz_mutant)
#				sys.stdout.write("\n")
#
#			#if dGcode[sz_codon][0] != dGcode[sz_mutant][0]:
#			#	sys.stdout.write("\t".join(ll))
#			#	sys.stdout.write("\t" + sz_codon + "\t" + ",".join(dGcode[sz_codon]))
#			#	sys.stdout.write("\t" + sz_mutant + "\t" + ",".join(dGcode[sz_mutant]))
#			#	sys.stdout.write("\n")
#		else:
#			sys.stdout.write("\t".join(ll))
#			#sys.stdout.write("\t" + sz_codon + "\t" + ",".join(dGcode[sz_codon]))
#			#sys.stdout.write("\t" + sz_mutant + "\t" + ",".join(dGcode[sz_mutant]))
#			sys.stdout.write("\t" + sz_codon + "\t" + ",".join(l_wtcode))
#			sys.stdout.write("\t" + "|".join(l_mutantCodon) + "\t" + sz_mutant)
#			sys.stdout.write("\n")
		
		
#	else:
#		ll.insert(8, str(float(n_hqhits)/float(n_hits)))
#		# simplify columns for exon, gene_id and transcript_id
#		ll[9] = "exon"
#		ltmp = list(set(ll[10].split(",")))
#		ll[10] = ",".join(ltmp)
#		ltmp = list(set(ll[11].split(",")))
#		ll[11] = ",".join(ltmp)
					
#		if (args.b_genes_only and b_inGene and b_print) or (args.b_printAll and b_print):
#		if args.b_printAll and b_print:
#		if b_print:
#			sys.stdout.write("\t".join(ll) + "\t")
#			sys.stdout.write("\t".join(["-", "-", "-", "-"]))
#			sys.stdout.write("\n")
#		ll += ["-", "-", "-", "-"]

	# append sample results
#	ll.append(",".join(sample_genotypes))
#	ll.append(",".join(map(str, sample_depths)))
#	ll.append(",".join(map(str, sample_altRatios)))

		
	if b_print:
		sys.stdout.write("\t".join(lout) + "\n")

ffa.close()
