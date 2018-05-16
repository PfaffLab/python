#!/usr/bin/env python
#==============================================================================
# ?.py
#
# shawn driscoll
# 20120719
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Hopefully this ends up being some type of junction classification script
#==============================================================================

import sys,argparse,re
import HTSeq as hts
import pybedtools as pbr
import numpy as np
import hashlib

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Classify detected junctions against a reference annotation.")
parser.add_argument('samples',type=str,help="sample file of junctions you wish to classify.")
parser.add_argument('reference',type=str,help="GTF reference file.")
parser.add_argument('-o',type=str,dest="prefix",default="proto",help="Prefix for output files (default: proto)")
parser.add_argument('-v',dest="b_verbose",action="store_const",const=True,default=False,help="Verbose output (default: false)")

args = parser.parse_args()

# check for file
try:
	fin = open(args.samples,"r")
except IOError as e:
	print "Error: unable to find/open input file: " + args.samples
	sys.exit()

try:
	fin = open(args.reference,"r")
except IOError as e:
	print "Error: unable to find/open input file: " + args.reference
	sys.exit()

#==============================================================================
# variables
#==============================================================================

dGtf = {}
dJuncs = {}
lHits = [False,False]
lHitAligned = [False,False]
lHitIndex = [-1,-1]
lSuffix = [".ja.txt",".ta.txt",".jaraw.txt"]
nri = 0

#==============================================================================
# main script
#==============================================================================


# parse the GTF into a dict to make looking up transcripts easy later on
# NOTE: i should verify the exons are sorted properly at this point

if args.b_verbose:
	sys.stderr.write("> parsing GTF file...\n")

gff = hts.GFF_Reader(args.reference)
nri = 0
for feature in gff:
	if feature.type == "exon":
		nri += 1
		if args.b_verbose:
			if nri % 2048 == 0:
				sys.stderr.write("\r> features parsed: %d  " % nri)

		szTid = feature.attr['transcript_id']
		if szTid not in dGtf:
			dGtf[szTid] = {}
			dGtf[szTid]['features'] = []
			dGtf[szTid]['junctions'] = []
			dGtf[szTid]['strand'] = feature.iv.strand
			dGtf[szTid]['gene_id'] = ""
			if "gene_id" in feature.attr:
				dGtf[szTid]['gene_id'] = feature.attr['gene_id']

		if len(dGtf[szTid]['features']) > 0:
			if feature.iv.start < dGtf[szTid]['features'][0].iv.start:
				dGtf[szTid]['features'].insert(0, feature)
			else:
				dGtf[szTid]['features'].append(feature)
		else:
			dGtf[szTid]['features'].append(feature)
				

		# dGtf[szTid]['features'].append(feature)
		# update junction map
		if len(dGtf[szTid]['features']) > 1:
			dGtf[szTid]['junctions'].append(0)

if args.b_verbose:
	sys.stderr.write("\r> features parsed %d\n" % nri)
	sys.stderr.write("> intersecting junctions with GTF file...\n")

# intersect the sample junctions with the GTF file via pybedtools
a = pbr.BedTool(args.samples)
b = pbr.BedTool(args.reference)
lres = a.intersect(b,wao=True)

if args.b_verbose:
	sys.stderr.write("> parsing intersections...\n")

# iterate through intersections and build a dict of the junctions
for ri in lres:
	if ri.fields[14] == "exon":
		szJid = ri.fields[3]

		if szJid not in dJuncs:
			# new junction encountered, record some info
			dJuncs[szJid] = {}
			#dJuncs[szJid]['feature'] = hts.GenomicInterval(ri.fields[0],int(ri.fields[1]),int(ri.fields[2]),".")
			dJuncs[szJid]['count'] = ri.fields[4]
			arr = map(int,ri.fields[10].split(","))
			dJuncs[szJid]['intron'] = hts.GenomicInterval(ri.fields[0], int(ri.fields[1])+arr[0], int(ri.fields[2])-arr[1], ".")
			dJuncs[szJid]['left'] = hts.GenomicInterval(ri.fields[0], int(ri.fields[1]), int(ri.fields[1])+arr[0], ".")
			dJuncs[szJid]['right'] = hts.GenomicInterval(ri.fields[0], int(ri.fields[2])-arr[1], int(ri.fields[2]), ".")
			dJuncs[szJid]["transcripts"] = []
			dJuncs[szJid]["transcript_results"] = {}

		m = re.search('transcript_id \"([\w\s\d\@\#\$\%\^\&\*\!\;\:\_\-\+\=\,\.\|\\\/\(\)\[\]]+)\"', ri.fields[20])
		# append transcript
		if m:
			szTid = m.group(1)
		else:
			szTid = "unknown"

		# don't add in entries for null hits - remember the intersection is set to include elements in the
		# junctions file that did not hit any annotation features.
		if szTid not in dJuncs[szJid]['transcripts'] and szTid != "unknown":
			dJuncs[szJid]['transcripts'].append(szTid)

		#dJuncs[szJid]['transcripts'][szTid].append(hts.GenomicFeature(szTid,"exon",hts.GenomicInterval(ri.fields[0],int(ri.fields[15])-1,int(ri.fields[16]),ri.fields[18])))

# now everything should be organized for analysis
# now we can loop through the junctions dict and evaluate each junction

#if args.b_verbose:
#	sys.stderr.write("> evaluating junctions: ")

nri = 0
for szJid in dJuncs.keys():
	# how many transcripts did this junction hit?
	lj = dJuncs[szJid]
	nTranscripts = len(lj['transcripts'])

	nri += 1
	if args.b_verbose:
		if nri % 128 == 0:
			sys.stderr.write("\r> evaluating junctions: %d" % nri)

	if nTranscripts > 0:
		# look at each transcript and evaluate the condition of the junction.
		for szTid in lj['transcripts']:

			# loop through the exon features of this transcript and find the hits for the left and
			# right sides of the junction

			lHits = [False,False]
			lHitIndex = [-1,-1]
			lLeftBetween = [-1,-1]
			lRightBetween = [-1,-1]
			lHitAligned = [False,False]

			# get transcript record
			lt = dGtf[szTid]
			nExons = len(lt['features'])
			nJunctions = nExons-1
			nAnnotatedJunctionIndex = -1
			#sys.stderr.write(str(nExons) + "\t" +  str(len(lt['junctions'])) + "\n")

			# start a dict entry for the results of this transcript in terms of the
			# current junction
			dJuncs[szJid]['transcript_results'][szTid] = {}

			for i in range(len(lt['features'])):
				# check for overlaps
				if lt['features'][i].iv.overlaps(lj['left']):
					lHits[0] = True
					lHitIndex[0] = i
					if lj['left'].end == lt['features'][i].iv.end:
						lHitAligned[0] = True

				if lt['features'][i].iv.overlaps(lj['right']):
					lHits[1] = True
					lHitIndex[1] = i
					if lj['right'].start == lt['features'][i].iv.start:
						lHitAligned[1] = True

				# check for upstream/downstream from entire isoform
				#if i == 0:
				#	if lj['left'].end < lt['features'][i].iv.start:
				#		lLeftBetween = [-1,0]
				#elif i == nExons-1:
				#	if lj['right'].start > lt['features'][i].iv.end:
				#		lRightBetween = [nExons-1,-1]

				# check between the current exon and the last one
				if i == 0:
					if lj['left'].end < lt['features'][i].iv.start:
						lLeftBetween = [-1,0]
				elif i == nExons-1:
					if lj['right'].start > lt['features'][i].iv.end:
						lRightBetween = [nExons-1,-1]
				else:
					if lj['left'].start > lt['features'][i-1].iv.end and lj['left'].end < lt['features'][i].iv.start:
						lLeftBetween = [i-1,i]
					if lj['right'].start > lt['features'][i-1].iv.end and lj['right'].end < lt['features'][i].iv.start:
						lRightBetween = [i-1,i]

			# save back to results
			dJuncs[szJid]['transcript_results'][szTid]['hits'] = [lHits[x] for x in range(2)]
			dJuncs[szJid]['transcript_results'][szTid]['hit_index'] = [lHitIndex[x] for x in range(2)]
			dJuncs[szJid]['transcript_results'][szTid]['left_between'] = [lLeftBetween[x] for x in range(2)]
			dJuncs[szJid]['transcript_results'][szTid]['right_between'] = [lRightBetween[x] for x in range(2)]
			dJuncs[szJid]['transcript_results'][szTid]['hit_aligned'] = [lHitAligned[x] for x in range(2)]
			dJuncs[szJid]['transcript_results'][szTid]['annotated'] = False
			dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = None
			dJuncs[szJid]['transcript_results'][szTid]['description'] = ""

			# check conditions of hits
			szType = ""
			if lHits[0] and lHits[1]:
				# both ends of the junction hit annotated exons - this could mean several things
				szType = "exon,exon;"
				if lHitIndex[1]-lHitIndex[0] == 1:
					# assign comparison junction index. since the exons hit are adjacent whether
					# the hits are aligned or not, this exon-exon junction is the one to compare
					# this junction to
					dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[0]

					if lHitAligned[0] and lHitAligned[1]:
						# this is a match
						if nExons > 1:
							try:
								dGtf[szTid]['junctions'][lHitIndex[0]] = lj['count']
								szType += "annotated;"
								dJuncs[szJid]['transcript_results'][szTid]['annotated'] = True
							except IndexError as e:
								sys.stderr.write("warning: index assignment out of range")
					else:
						# one or both sides are misaligned
						if not lHitAligned[0]:
							szType += "lma;"
						if not lHitAligned[1]:
							szType += "rma;"

				elif lHitIndex[1] - lHitIndex[0] > 1:
					# this is exon to exon with exon exclusion
					szType += "exon_exclusion;"

					# if the left or right side of this junction is misaligned it changes
					# the junction we want to comapre it to in the reference.
					if not lHitAligned[0] and lHitAligned[1]:
						# right side is aligned we want junction at that exon
						dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[1]-1
					elif lHitAligned[0] and not lHitAligned[1]:
						# left side is aligned so we want the junction that stems (or terminates) at that exon
						dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[0]
					elif (lHitAligned[0] and lHitAligned[1]) or (not lHitAligned[0] and not lHitAligned[1]):
						# if we have the same condition at both ends then we'll take strand into account. for 
						# positive strand we'll take the junction that starts on the left side and for
						# negative strand we'll take the junction that starts from the right
						if lt['features'][i].iv.strand == "+":
							dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[0]
						else:
							dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[1]-1
					else:
						sys.stderr.write("warning: unexpected condition @ exon to exon with exon exclusion\n")

					# one or both sides are misaligned
					if not lHitAligned[0]:
						szType += "lma;"
					if not lHitAligned[1]:
						szType += "rma;"

				elif lHitIndex[1] == lHitIndex[0]:
					# this is a junction completely contained within an exon (sometimes this happens with
					# 3' or 5' exons)
					szType += "intra_exon;"

			elif not lHits[0] and lHits[1]:
				# evaluate conditions where the right side hit but the left did not
				szType = "intron,exon;"

				# if the right side is hitting exon 0 then there is not a junction
				# to compare this one to. since this value is set to None by default
				# we only need to set it if the exon hit is any other exon
				if dJuncs[szJid]['transcript_results'][szTid]['hit_index'][1] > 0:
					# comparison junction will be the one that hits the same exon as the right side-1
					dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[1]-1
				
				# is the right side aligned?
				if lHitAligned[1]:
					szType += "aligned;"
				else:
					szType += "misaligned;"

				# where does the left side fall
				# left side is between two exons of the isoform. does this exlude any exons?
				if lLeftBetween[1] < lHitIndex[1]:
					szType += "exon_exclusion;"
				else:
					szType += "adjacent_novel_exon;"

				# does the junction extend past the end of this isoform?
				if lLeftBetween[0] == -1 and lHitIndex[1] != 0:
					if lt['strand'] == "+":
						szType += "5'_lost;"
					else:
						szType += "3'_lost;"


			elif lHits[0] and not lHits[1]:
				# evaluate conditions where the left side hit but the right did not
				szType = "exon,intron;"

				# if the left side hits the last exon in the transcript then there
				# is no junction to compare it to. this value is set to None by default
				# so we only need to set it if the exon index is at any other index.
				if dJuncs[szJid]['transcript_results'][szTid]['hit_index'][0] < (len(dGtf[szTid]['features'])-1):
					# comparison junction will be the one that hits the same exon as the left side
					dJuncs[szJid]['transcript_results'][szTid]['comp_index'] = lHitIndex[0]
				
				# is the right side aligned?
				if lHitAligned[0]:
					szType += "aligned;"
				else:
					szType += "misaligned;"

				# where does right side fall?
				if lRightBetween[0] > lHitIndex[0]:
					szType += "exon_exclusion;"
				else:
					szType += "adjacent_novel_exon;"

				# does the junction extend past the end of this isoform?
				if lRightBetween[1] == -1 and lHitIndex[0] != nExons-1:
					if lt['strand'] == "+":
						szType += "3'_lost;"
					else:
						szType += "5'_lost;"

			elif not lHits[0] and not lHits[1]:
				szType = "intron,intron;"

			else:
				sys.stderr.write("warning: unexpected landing condition\n")

			# comparison junction will be the one that hits the same exon as the right side-1
			dJuncs[szJid]['transcript_results'][szTid]['description'] = szType

			#sys.stdout.write(szJid + "\t")
			#sys.stdout.write(szTid + "\t")
			#sys.stdout.write(",".join(map(str,lHitIndex)) + "\t")
			#sys.stdout.write(",".join(map(str,lHitAligned)) + "\t")
			#sys.stdout.write(",".join(map(str,lLeftBetween)) + "\t")
			#sys.stdout.write(",".join(map(str,lRightBetween)))
			#sys.stdout.write("\t" + szType)
			#sys.stdout.write("\t" + szType + "\t" + str(nExons) + "\t" + str(lj['count']))
			#sys.stdout.write("\n")

# now we have to loop back through the junctions again since we have quantified
# all of the junctions for the annotation. now we can print out junctions and their 
# compairson junction as well as attempt to select the one, of many possible, transcript
# this that's the best match for the junction
if args.b_verbose:
	sys.stderr.write("\r> evaluating junctions: %d\n" % nri)
	sys.stderr.write("> writing junction analysis results...\n")

# output is written to the first file
fout = open(args.prefix + lSuffix[0],"w")

for szJid in sorted(dJuncs.keys()):
	# how many transcripts did this entry hit?
	lj = dJuncs[szJid]
	nTranscripts = len(lj['transcript_results'])
	# make a string location out of the intron 
	szLoc = lj['intron'].chrom + ":" + str(lj['intron'].start) + "-" + str(lj['intron'].end)

	if nTranscripts > 0:

		lMatches = []

		# first off is there a match?
		for szTid in lj['transcript_results'].keys():
			if lj['transcript_results'][szTid]['annotated']:
				lMatches.append(szTid)

		nCount = lj['count']
		nCompCount = 0
		nCountRatio = 0
		if len(lMatches) > 0:
			# this junction matches an annotated junction for at least one of the
			# isoforms hit, so it's done
			nCompCount = dGtf[lMatches[0]]['junctions'][lj['transcript_results'][lMatches[0]]['comp_index']]
			nCountRatio = 1

			fout.write(szJid)
			fout.write("\t" + szLoc)
			fout.write("\t" + hashlib.md5(szLoc).hexdigest())
			fout.write("\t" + ",".join(lMatches))
			fout.write("\t=")
			fout.write("\t" + str(lj['count']))
			fout.write("\t" + str(nCompCount))
			fout.write("\t" + str(nCountRatio))
			fout.write("\t" + lj['transcript_results'][lMatches[0]]['description'])
			fout.write("\n")
		else:
			# novel junction because there are no matches - now figure out which transcript it is the best match
			# to

			# using a simple score system we should be able to figure out which isoform is the best match
			# best: exon to exon, aligned
			# for a missing exon we subtract 2, for a missing alignment we subtract 1
			# so a exon to intron with an aligned junction at the exon has a score of 2
			# but with no alignment at the exon we get 1.
			arScores = np.array([4 for x in range(nTranscripts)])
			# loop through transcripts
			ni = 0
			lKeys = lj['transcript_results'].keys()
			for szTid in lKeys:
				if not lj['transcript_results'][szTid]['hits'][0]:
					arScores[ni] -= 2
				elif not lj['transcript_results'][szTid]['hit_aligned'][0]:
					arScores[ni] -= 1

				if not lj['transcript_results'][szTid]['hits'][1]:
					arScores[ni] -= 2
				elif not lj['transcript_results'][szTid]['hit_aligned'][1]:
					arScores[ni] -= 1
				ni += 1

			nMaxPos = arScores.argmax()
			nMaxScore = arScores.max()

			#if szJid == "JUNC00000412":
			#	print lKeys,arScores
			#	sys.exit()

			szTid = lKeys[nMaxPos]

			if lj['transcript_results'][szTid]['comp_index'] is None:
				nCompCount = "na"
				nCountRatio = "inf"
			else:
				try:
					nCompCount = dGtf[szTid]['junctions'][lj['transcript_results'][szTid]['comp_index']]
					if nCompCount == 0:
						nCountRatio = "inf"
					else:
						nCountRatio = float(nCount)/float(nCompCount)
				except:
					nCompCount = -1
					nCountRatio = "na"

			fout.write(szJid)
			fout.write("\t" + szLoc)
			fout.write("\t" + hashlib.md5(szLoc).hexdigest())			
			fout.write("\t" + szTid)
			fout.write("\tnovel")
			fout.write("\t" + str(lj['count']))
			fout.write("\t" + str(nCompCount))
			fout.write("\t" + str(nCountRatio))
			fout.write("\t" + lj['transcript_results'][szTid]['description'])
			fout.write("\n")

			#sys.stdout.write(szJid)
			#sys.stdout.write("\t" + " ")
			#sys.stdout.write("\t" + str(lj['count']))
			#sys.stdout.write("\t" + str(dGtf[szTid]['junctions'][lj['transcript_results'][szTid]['comp_index']]))
			#sys.stdout.write("\n")

fout.close()

if args.b_verbose:
	sys.stderr.write("> writing transcript results...\n")

# write transcript results out to second file
fout = open(args.prefix + lSuffix[1],"w")

szCov = "="

# print out all of the isoform junction count summaries
for szTid in sorted(dGtf.keys()):

	# is this transcript's junction set totally covered?
	szCov = "="
	if len(dGtf[szTid]['features']) > 1:
		if sum(map(int,dGtf[szTid]['junctions'])) == 0:
			szCov = "n"
		elif min(map(int,dGtf[szTid]['junctions'])) == 0:
			szCov = "i"
	else:
		szCov = "nj"

	fout.write(szTid)
	fout.write("\t" + dGtf[szTid]['strand'])
	fout.write("\t" + szCov)
	if len(dGtf[szTid]['features']) > 1:
		fout.write("\t" + ",".join(map(str,dGtf[szTid]['junctions'])))
	else:
		fout.write("\tna")

	fout.write("\n")

fout.close()


if args.b_verbose:
	sys.stderr.write("> done!\n")




	#print szJid + "\t" + str(nTranscripts) + "\t" + ",".join(dJuncs[szJid]['transcripts'].keys())


