#!/usr/bin/env python
#==============================================================================
# coverage-analysis.py
#
# Shawn Driscoll
# 20120705
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Parses output of coverageBed against a GTF file and produces a report
# of the per-exon coverage and hit counts as well as exon "expressions".
# input via stdin must be sorted on id column and by position 
# sort -k9,9 -k4,4n will do it - also make sure that only exon rows are 
# passed to this script.
#==============================================================================

import sys,re,math

#==============================================================================
# variables
#==============================================================================

# create a transcript's dictionary
dTranscripts = {}
szTid = ""
szGid = ""
szLastTid = ""
n_hitThreshold = 0
sz_class = ""

#==============================================================================
# functions
#==============================================================================

#
# sort
# Sorts values in da and returns a vector of the new order of indices
# as a result of the sorting
def sort(da):
	# length 
	nlen = len(da)
	# create indices
	dindex = range(nlen)
	# flags
	bUnsorted = True
	ni = 0
	ntmp = 0

	while bUnsorted:
		bUnsorted = False
		ni = 0
		while not bUnsorted and ni < nlen-1:
			if da[ni] > da[ni+1]:
				# rearrange values
				ntmp = da[ni]
				da[ni] = da[ni+1]
				da[ni+1] = ntmp
				# rearrange indices
				ntmp = dindex[ni]
				dindex[ni] = dindex[ni+1]
				dindex[ni+1] = ntmp
				# set flag
				bUnsorted = True

			ni += 1

	return dindex


#==============================================================================
# main script
#==============================================================================

# read from stdin
sys.stderr.write("parsing coveragebed output on stdin...\n")
fin = sys.stdin

for szl in fin:
	szl = szl.strip()
	arl = szl.split("\t")

	# parse out transcript id and gene id
	#m = re.search('gene_id \"([\w\s\d\@\#\$\%\^\&\*\!\:\;\_\-\+\=\,\.\|\\\/\(\)\[\]]+)\"', arl[8])
	m = re.search('gene_id \"([^\"]+)\"', arl[8])
	szGid = ""
	if m:
		szGid = m.group(1)

	m = re.search('transcript_id \"([^\"]+)\"', arl[8])
	szTid = ""
	if m:
		szTid = m.group(1)

	if szTid not in dTranscripts:
		# create new entry in dictionary
		dTranscripts[szTid] = {'ecount':0, 'ehit': 0, 'ehits':[], 'ecov':[], 'elen':[], 'estarts':[], 'strand':arl[6], 'name':szGid, 'hits':0, 'length':0, 'effectiveLength':0}

	# insert information into dictionary
	dTranscripts[szTid]['ecount'] += 1
	dTranscripts[szTid]['ehits'].append(arl[9])
	if int(arl[9]) > n_hitThreshold:
		dTranscripts[szTid]['ehit'] += 1

	dTranscripts[szTid]['ecov'].append(arl[12])
	dTranscripts[szTid]['elen'].append(arl[11])
	dTranscripts[szTid]['estarts'].append(int(arl[3]))
	dTranscripts[szTid]['hits'] += int(arl[9])
	dTranscripts[szTid]['length'] += int(arl[11])
	dTranscripts[szTid]['effectiveLength'] += int(arl[11]) * float(arl[12])

	szLastTid = szTid

# make sure exons are sorted properly
for szTid in sorted(dTranscripts.keys()):
	# sort indicies by exon start positions then apply index order to 
	# each of the value lists for the transcript
	dsi = sort(dTranscripts[szTid]['estarts'])
	nlen = len(dsi)
	dStarts = range(nlen)
	dHits = range(nlen)
	dCov = range(nlen)
	dLen = range(nlen)

	# reorder values
	for i in range(nlen):
		dHits[i] = dTranscripts[szTid]['ehits'][dsi[i]]
		dStarts[i] = dTranscripts[szTid]['estarts'][dsi[i]]
		dCov[i] = dTranscripts[szTid]['ecov'][dsi[i]]
		dLen[i] = dTranscripts[szTid]['elen'][dsi[i]]

	# copy back into vector
	for i in range(nlen):
		dTranscripts[szTid]['ehits'][i] = dHits[i]
		dTranscripts[szTid]['estarts'][i] = dStarts[i]
		dTranscripts[szTid]['ecov'][i] = dCov[i]
		dTranscripts[szTid]['elen'][i] = dLen[i]


# print it all out

# print header
print "\t".join(['transcript_id','gene_id','strand','length', 'effectiveLength','hits','num_exons','num_hit_exons','exon_lengths','exon_coverages','exon_hits','exon_starts'])

# print transcripts
for szTid in sorted(dTranscripts.keys()):
	sz_class = "="
	if min(map(int,dTranscripts[szTid]['ehits'])) < n_hitThreshold:
		sz_class = "i"

	sys.stdout.write(szTid + "\t" + dTranscripts[szTid]['name'] + "\t" + dTranscripts[szTid]['strand'] + "\t")
	sys.stdout.write(str(dTranscripts[szTid]['length']) + "\t" + str(math.floor(dTranscripts[szTid]['effectiveLength'])) + "\t" + str(dTranscripts[szTid]['hits']) + "\t")
	sys.stdout.write(sz_class + "\t")
	sys.stdout.write(str(dTranscripts[szTid]['ecount']) + "\t" + str(int(math.floor(dTranscripts[szTid]['ehit']))) + "\t")
	sys.stdout.write(",".join(dTranscripts[szTid]['elen']) + "\t")
	sys.stdout.write(",".join(dTranscripts[szTid]['ecov']) + "\t")
	sys.stdout.write(",".join(dTranscripts[szTid]['ehits']) + "\t")
	sys.stdout.write(",".join(map(str,dTranscripts[szTid]['estarts'])) + "\n")


