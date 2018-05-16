#!/usr/bin/env python
#==============================================================================
# vcf-get-genes.py
#
# Shawn Driscoll
# 20120430
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Fetches gene names for the SNPS in the passed vcf file
#==============================================================================

import sys,argparse,hashlib,re
import pybedtools as bpr

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser("Import gene and transcript names from a GTF file into a VCF file.")

parser.add_argument('vcf', action="store",type=str,help="VCF format SNP file")
parser.add_argument("gtf", action="store", type=str, help="GTF gene annotation")
parser.add_argument("-vcf", dest="outVcf", action="store_const", const=True, default=False, help="Output VCF format (default: parsed output)")

args = parser.parse_args()

#==============================================================================
# setup variables
#==============================================================================

dSnps = {}
dGids = {}
dTids = {}


#==============================================================================
# main function
#==============================================================================


# parse the VCF file into memory so we have a record of what all is there
fin = open(args.vcf,"r")

for szl in fin:
	szl = szl.strip()
	arl = szl.split("\t")

	# create hash for this SNP based on location
	szh = hashlib.md5(arl[0] + arl[1]).hexdigest()

	if szh not in dSnps:
		dSnps[szh] = arl
	else:
		sys.stderr.write("> duplicate SNP at " + arl[0] + ":" + arl[1] + "\n")

fin.close()

# intersect VCF with the GTF and parse the results

a = bpr.BedTool(args.vcf)
b = bpr.BedTool(args.gtf)

res = a.intersect(b,wao=True)

# iterate through intersection. we'll have 11 columns of VCF and the rest are GTF

for ri in res:
	szh = hashlib.md5(ri.fields[0] + ri.fields[1]).hexdigest()

	if ri.fields[13] == "exon":
		# find an intersection, parse out the gene name and transcript name

		szGid = ""
		szTid = ""

		m = re.search('gene_id \"([^\"]+)\"', ri.fields[19])
		if m:
			szGid = m.group(1)
		else:
			sys.stderr.write("> error: no gene found: (" + szl + ")\n")

		# parse out transcript id
		m = re.search('transcript_id \"([\w\s\d\@\#\$\%\^\&\*\!\;\:\_\-\+\=\,\.\|\\\/\(\)\[\]]+)\"', ri.fields[19])
		if m:
			szTid = m.group(1)
		else:
			sys.stderr.write("> error: no transcript id found: (" + szl + ")\n")

		# append gene
		if szh not in dGids:
			dGids[szh] = []

		dGids[szh].append(szGid)

		# append transcript
		if szh not in dTids:
			dTids[szh] = []

		dTids[szh].append(szTid)


# intersection is complete now we can print the original list back out with any 
# discovered gene info

for szh in dSnps.keys():
	
	arSnp = dSnps[szh]

	szGid = ""
	szTid = ""

	if szh in dGids:
		arGids = list(set(dGids[szh]))
		arTids = list(set(dTids[szh]))

	if args.outVcf:
		# append gene name
		arSnp[7] = arSnp[7] + ";gene_id=" + ",".join(arGids) + ";transcript_id=" + ",".join(arTids)
		print "\t".join(arSnp)
	else:
		arOut = arSnp[0:7]
		# parse slot 7 for the DP and DP4 arguments
		sztmp = ""
		
		m = re.search('DP\=([0-9]+)\;',arSnp[7])
		if m:
			sztmp = m.group(1)

		arOut.append(sztmp)

		sztmp = ""
		m = re.search('DP4\=([0-9\,]+)\;',arSnp[7])
		if m:
			sztmp = m.group(1)

		# split list - first two are reference alleal and second two are mutant alleal
		artmp = map(int,sztmp.split(","))
		arOut.append(str(artmp[0]+artmp[1]))
		arOut.append(str(artmp[2]+artmp[3]))

		print "\t".join(arOut) + "\t" + ",".join(arGids) + "\t" + ",".join(arTids)

# done



