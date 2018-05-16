#!/usr/bin/env python
#==============================================================================
# annot-peaks.py
#
# Shawn Driscoll
# 20150203
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Script to annotate chip peaks as exported from Homer (or optionally a bed file)
#==============================================================================

import sys, argparse, re

# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# main
#==============================================================================

_HBIN = 16000

def main(args):

	# variables

	# parse GTF into lookup hashes
	(exons, tss, introns, tinfo, cds, pro) = load_gtf(args.gtf_file, False)

	# open peaks
	try:
		fin = open(args.peaks_file, "r")
	except IOError as e:
		sys.stderr.write("IOError ({}): {}\n".format(e.errno, e.strerror))
		return 1

	if args.bed:
		# peak into file
		szl = fin.readline()
		aln = szl.strip().split("\t")
		nbedcols = len(aln)
		# print header out
		haln = ["chrom", "start", "end"]
		if nbedcols > 3:
			haln += ["name"]
		if nbedcols > 4:
			haln += ["score"]
		# that's enough
		nbedcols = 5

		haln += ["withinFeature", "withinGene", "withinTranscript", "nearestDsGene", "nearestDsTranscript", "DsDistance", 
				"nearestUsGene", "nearestUsTranscript", "UsDistance"]

		print "\t".join(haln)
		fin.close()
		fin = open(args.peaks_file, "r")

		# set some indices
		idx_chrom = 0
		idx_pstart = 1
		idx_pend = 2

	else:

		idx_chrom = 1
		idx_pstart = 2
		idx_pend = 3


	for szl in fin:

		if not args.bed:
			if re.search("^\#PeakID", szl) or re.search("Start", szl):
				# column headers, edit this and print it back out
				# get transcript information fields
	#			ttmp = tinfo[tinfo.keys()[0]]
	#			ttmp_names = sorted(ttmp.keys())
				aln = szl.strip().split("\t")
				aln += ["withinFeature", "withinGene", "withinTranscript", "nearestDsGene", "nearestDsTranscript", "DsDistance", 
					"nearestUsGene", "nearestUsTranscript", "UsDistance"]

				print "\t".join(aln)

				# hinfo += [dgn, dbest_id, dbest, ugn, ubest_id, ubest]
				continue

			elif szl[0] == "#":
				continue

		aln = szl.strip().split("\t")

		# peak range
		pkr = [int(aln[idx_pstart]), int(aln[idx_pend])]
		if pkr[1]-pkr[0] > 1:
			pkcent = round(float(pkr[0]+pkr[1])/2)
		else:
			pkcent = pkr[1]
	
		eres = lookup_pos(exons, aln[idx_chrom], pkcent)
		ires = lookup_pos(introns, aln[idx_chrom], pkcent)
		pres = lookup_pos(pro, aln[idx_chrom], pkcent)

		phit = []
		if len(pres) > 0:
			for i in range(len(pres)):
				if overlap(pres[i][1:3], [pkcent, pkcent])==1:
					phit.append(list(pres[i]))

		ehit = []
		if len(eres) > 0:
			for i in range(len(eres)):
				if overlap(eres[i][1:3], [pkcent, pkcent])==1:
					ehit.append(list(eres[i]))

		ihit = []
		if len(ires) > 0:
			for i in range(len(ires)):
				if overlap(ires[i][1:3], [pkcent, pkcent])==1:
					ihit.append(list(ires[i]))


		if len(phit) > 0:
			# make info
			tid = {}
			gn = {}
			for i in range(len(phit)):
				ttmp = phit[i][0].split(",")
				for t in ttmp:
					tid[t] = 0
					gn[tinfo[t]['gene_name']] = 0

			hinfo = ["promoter", ",".join(gn.keys()), ",".join(tid.keys())]

		elif len(ehit) > 0:
			# make info
			tid = {}
			gn = {}
			for i in range(len(ehit)):
				ttmp = ehit[i][0].split(",")
				for t in ttmp:
					tid[t] = 0
					gn[tinfo[t]['gene_name']] = 0

			hinfo = ["exon", ",".join(gn.keys()), ",".join(tid.keys())]

		elif len(ihit) > 0:
			# make info
			tid = {}
			gn = {}
			for i in range(len(ihit)):
				ttmp = ihit[i][0].split(",")
				for t in ttmp:
					tid[t] = 0
					gn[tinfo[t]['gene_name']] = 0

			hinfo = ["intron", ",".join(gn.keys()), ",".join(tid.keys())]
		else:
			hinfo = ["intergenic", "-", "-"]

		# find nearest tss

		k0 = make_hbin(aln[idx_chrom], pkcent, pkcent)[0]

		dbest = 1e9
		dbest_id = ""
		ubest = -1e9
		ubest_id = ""

		if k0 in tss:
			# find nearest in both directions
			dbest = 1e9
			dbest_id = ""
			ubest = -1e9
			ubest_id = ""

			for i in range(len(tss[k0])):
				diff = int(tss[k0][i][1])-pkcent
				if diff < 0 and diff > ubest:
					ubest = diff
					ubest_id = tss[k0][i][0]
				elif diff > 0 and diff < dbest:
					dbest = diff
					dbest_id = tss[k0][i][0]

		if dbest_id == "":
			# find nearest downstream tss
			for j in range(1,int(round(1e6/_HBIN))):
				k = make_hbin(aln[idx_chrom], pkcent+j*_HBIN, pkcent+j*_HBIN)[0]
				if k in tss:
					# got it
					for i in range(len(tss[k])):
						diff = int(tss[k][i][1])-pkcent
						if diff < dbest:
							dbest = diff
							dbest_id = tss[k][i][0]

					break

		if ubest_id == "":
			# find nearest upstream tss
			for j in range(1,int(round(1e6/_HBIN))):
				if pkcent - j*_HBIN < 0:
					break

				k = make_hbin(aln[idx_chrom], pkcent-j*_HBIN, pkcent-j*_HBIN)[0]
				if k in tss:
					# got it
					for i in range(len(tss[k])):
						diff = int(tss[k][i][1])-pkcent
						if diff > dbest:
							dbest = diff
							dbest_id = tss[k][i][0]

					break

		if dbest_id != "":
			ttmp = dbest_id.split(",")
			gn = {}
			for t in ttmp:
				gn[tinfo[t]['gene_name']] = 0

			dgn = ",".join(gn.keys())
			dbest = int(dbest)
		else:
			dgn = '-'
			dbest_id = '-'
			dbest = '-'

		if ubest_id != "":
			ttmp = ubest_id.split(",")
			gn = {}
			for t in ttmp:
				gn[tinfo[t]['gene_name']] = 0

			ugn = ",".join(gn.keys())
			ubest = int(ubest)
		else:
			ugn = '-'
			ubest_id = "-"
			ubest = "-"

		hinfo += [dgn, dbest_id, dbest, ugn, ubest_id, ubest]

		print "\t".join(aln + map(str, hinfo))

	fin.close()

	return 0


# --
# function to parse the GTF and create hashes of all features
def load_gtf(f, cds_tss_only):

	hexons = {}
	htss = {}
	hintrons = {}
	hcds = {} # hash of transcript names that have CDS regions in the annotation
	hpro = {} # promoters will be 1000 upstream / 100 downstream of tss

	nexons = 0
	ntss = 0
	nintrons = 0
	npro = 0

	tid = ""
	tid_last =""

	tbuffer = [] # buffer for rows of a single transcript as it is read
	tintrons = []
	attr = {}
	attr_last = {}

	tinfo = {}

	try:
		fin = open(f, "r")
	except IOError as e:
		sys.stderr.write("IO Error ({}): {}\n".format(e.errno, e.strerror))
		raise

	sys.stderr.write("[load_gtf] reading {}\n".format(f))

	# parse the file
	for szl in fin:
		aln = szl.strip().split("\t")

		# parse attribute field
		attr = parse_gtf_attr(aln[8])
		# add strand into attributes
		attr['strand'] = aln[6]
		if 'gene_name' not in attr:
			attr['gene_name'] = attr['transcript_id']

		if aln[2] == "CDS":
			hcds[attr['transcript_id']] = 0
		elif aln[2] != "exon":
			continue

		tid = attr['transcript_id']

		# add transcript information into hash
		if tid not in tinfo:
			tinfo[tid] = dict(attr)
		
		# new transcript reached so process the last one
		if tid != tid_last and tid_last != "":

			# deal with last
			
			if len(tbuffer) > 0:
				tintrons = []
				tstrand = tbuffer[0][6]

				# sort the buffer, add introns
				if len(tbuffer) > 0:
					tintrons = []
					tstrand = tbuffer[0][6]

					# sort the buffer, add introns
					if len(tbuffer) > 1:
						tbuffer.sort(key=lambda x:int(x[3]))
						# make introns
						for i in range(1, len(tbuffer)):
							# make list of tid and start/end of intron
							itmp = [attr_last['transcript_id'], int(tbuffer[i-1][4])+1, int(tbuffer[i][3])-1]
							# make hash bin key
							k0 = make_hbin(tbuffer[i][0], int(tbuffer[i-1][4])+1, int(tbuffer[i][3])-1)
							for k in k0:
								# add bin and intron to the introns hash
								idx = -1
								if k not in hintrons:
									hintrons[k] = []
								else:
									# bucket exists, see if this intron exists
									for j in range(len(hintrons[k])):
										if intron_cmp(hintrons[k][j], itmp) == 1:
											idx = j
											break

								# add to hash
								if idx >= 0:
									# append transcript id to the intron's transcript id string
									hintrons[k][idx][0] = hintrons[k][idx][0] + "," + itmp[0]
								else:
									nintrons += 1
									hintrons[k].append(list(itmp))

					# add exons
					for i in range(len(tbuffer)):
						etmp = [attr_last['transcript_id'], int(tbuffer[i][3]), int(tbuffer[i][4])]
						k0 = make_hbin(tbuffer[i][0], tbuffer[i][3], tbuffer[i][4])
						for k in k0:
							idx = -1
							if k not in hexons:
								hexons[k] = []
							else:
								# bucket exists, see if exon exists
								for j in range(len(hexons[k])):
									if intron_cmp(hexons[k][j], etmp)==1:
										# exon exists, append tid
										idx = j
										break

							# add to hash
							if idx >= 0:
								hexons[k][idx][0] = hexons[k][idx][0] + "," + etmp[0]
							else:
								nexons += 1
								hexons[k].append(list(etmp))

					# add tss to htss. if strand is - then this is tbuffer[-1][4]
					# else it is tbuffer[0][3]

					if tstrand == "+" or len(tbuffer)==1:
						ttmp = [attr_last['transcript_id'], int(tbuffer[0][3]), tstrand]
					elif tstrand == "-":
						ttmp = [attr_last['transcript_id'], int(tbuffer[-1][4]), tstrand]

					add_tss = True
					if cds_tss_only:
						if attr_last['transcript_id'] not in hcds:
							add_tss = False

					if add_tss:
						k0 = make_hbin(tbuffer[0][0], ttmp[1], ttmp[1])

						for k in k0:
							idx = -1
							if k not in htss:
								htss[k] = []
							else:
								# bucket exists, see if the tss is there already
								for j in range(len(htss[k])):
									if tss_cmp(htss[k][j], ttmp):
										idx = j
										break
							# add to hash
							if idx >= 0:
								# append transcript id
								htss[k][idx][0] = htss[k][idx][0] + "," + ttmp[0]
							else:
								ntss += 1
								htss[k].append(list(ttmp))


						# add promoter
						ptmp = [ttmp[0], ttmp[1], ttmp[1], tstrand]
						if tstrand == "+":
							ptmp[1] -= 1000
							ptmp[2] += 100
						else:
							ptmp[1] -= 100
							ptmp[2] += 1000

						k0 = make_hbin(tbuffer[0][0], ptmp[1], ptmp[2])

						for k in k0:
							idx = -1
							if k not in hpro:
								hpro[k] = []
							else:
								# bucket exists, see if the tss is there already
								for j in range(len(hpro[k])):
									if tss_cmp(hpro[k][j], ptmp):
										idx = j
										break
							# add to hash
							if idx >= 0:
								# append transcript id
								hpro[k][idx][0] = hpro[k][idx][0] + "," + ptmp[0]
							else:
								npro += 1
								hpro[k].append(list(ptmp))


			# clear the buffer for the next transcript
			tbuffer = []

		tbuffer.append(aln)
		tid_last = tid
		attr_last = dict(attr)

	fin.close()


	if len(tbuffer) > 0:
		tintrons = []
		tstrand = tbuffer[0][6]

		# sort the buffer, add introns
		if len(tbuffer) > 1:
			tbuffer.sort(key=lambda x:int(x[3]))
			# make introns
			for i in range(1, len(tbuffer)):
				# make list of tid and start/end of intron
				itmp = [attr_last['transcript_id'], int(tbuffer[i-1][4])+1, int(tbuffer[i][3])-1]
				# make hash bin key
				k0 = make_hbin(tbuffer[i][0], int(tbuffer[i-1][4])+1, int(tbuffer[i][3])-1)
				for k in k0:
					# add bin and intron to the introns hash
					idx = -1
					if k not in hintrons:
						hintrons[k] = []
					else:
						# bucket exists, see if this intron exists
						for j in range(len(hintrons[k])):
							if intron_cmp(hintrons[k][j], itmp) == 1:
								idx = j
								break

					# add to hash
					if idx >= 0:
						# append transcript id to the intron's transcript id string
						hintrons[k][idx][0] = hintrons[k][idx][0] + "," + itmp[0]
					else:
						nintrons += 1
						hintrons[k].append(list(itmp))

		# add exons
		for i in range(len(tbuffer)):
			etmp = [attr_last['transcript_id'], int(tbuffer[i][3]), int(tbuffer[i][4])]
			k0 = make_hbin(tbuffer[i][0], tbuffer[i][3], tbuffer[i][4])
			for k in k0:
				idx = -1
				if k not in hexons:
					hexons[k] = []
				else:
					# bucket exists, see if exon exists
					for j in range(len(hexons[k])):
						if intron_cmp(hexons[k][j], etmp)==1:
							# exon exists, append tid
							idx = j
							break

				# add to hash
				if idx >= 0:
					hexons[k][idx][0] = hexons[k][idx][0] + "," + etmp[0]
				else:
					nexons += 1
					hexons[k].append(list(etmp))

		# add tss to htss. if strand is - then this is tbuffer[-1][4]
		# else it is tbuffer[0][3]

		if tstrand == "+" or len(tbuffer)==1:
			ttmp = [attr_last['transcript_id'], int(tbuffer[0][3]), tstrand]
		elif tstrand == "-":
			ttmp = [attr_last['transcript_id'], int(tbuffer[-1][4]), tstrand]


		add_tss = False
		if cds_tss_only:
			if attr_last['transcript_id'] in hcds:
				add_tss = True
		else:
			add_tss = True

		if add_tss:
			k0 = make_hbin(tbuffer[0][0], ttmp[1], ttmp[1])
			for k in k0:
				idx = -1
				if k not in htss:
					htss[k] = []
				else:
					# bucket exists, see if the tss is there already
					for j in range(len(htss[k])):
						if tss_cmp(htss[k][j], ttmp):
							idx = j
							break
				# add to hash
				if idx >= 0:
					# append transcript id
					htss[k][idx][0] = htss[k][idx][0] + "," + ttmp[0]
				else:
					ntss += 1
					htss[k].append(list(ttmp))

			# add promoter
			ptmp = [ttmp[0], ttmp[1], ttmp[1], tstrand]
			if tstrand == "+":
				ptmp[1] -= 1000
				ptmp[2] += 100
			else:
				ptmp[1] -= 100
				ptmp[2] += 1000

			k0 = make_hbin(tbuffer[0][0], ptmp[1], ptmp[2])

			for k in k0:
				idx = -1
				if k not in hpro:
					hpro[k] = []
				else:
					# bucket exists, see if the tss is there already
					for j in range(len(hpro[k])):
						if tss_cmp(hpro[k][j], ptmp):
							idx = j
							break
				# add to hash
				if idx >= 0:
					# append transcript id
					hpro[k][idx][0] = hpro[k][idx][0] + "," + ptmp[0]
				else:
					npro += 1
					hpro[k].append(list(ptmp))


	sys.stderr.write("[load_gtf] parsed {} exons; {} introns; {} tss; {} promoters\n".format(nexons, nintrons, ntss, npro))

	return (hexons, htss, hintrons, tinfo, hcds, hpro)

# --
# return list of items in hash bucket from a lookup hash. ref and pos are hashed and 
# looked up in hash, h. if nothing is found then an empty list is returned
def lookup_pos(h, ref, pos):
	# make key
	k = make_hbin(ref, pos, pos)
	k = k[0]

	if k in h:
		return(h[k])

	return []

def intron_cmp(i1, i2):
	if i1[1]==i2[1] and i1[2]==i2[2]:
		return 1
	return 0

def overlap(r1, r2):
	r1 = map(float, r1)
	r2 = map(float, r2)

	if r1[0] <= r2[1] and r1[1] >= r2[0]:
		return 1
	return 0

def tss_cmp(t1, t2):
	if t1[1]==t2[1] and t1[2]==t2[2]:
		return 1
	return 0


# --
# create hash index from reference and position for lookup tables
def make_hbin(ref, pos0, pos1):
	
	pos = map(int, [pos0, pos1])

	bins = [p/_HBIN for p in pos]
	k = []

	for b in range(bins[0], bins[1]+1):
		k.append("{}:{}".format(ref, b))

	return(k)

# --
# parse the attributes field of a GTF row into a hash
def parse_gtf_attr(sz):
	# split string on ;
	s1 = sz.split(";")
	attr = {}

	for i in range(len(s1)):
		s1[i] = s1[i].strip()
		if len(s1[i]) > 0:
			m = re.search("^([^\s]+) \"([^\"]+)\"", s1[i])
			attr[m.group(1)] = m.group(2)

	return attr
		




#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="...")
parser.add_argument('gtf_file', type=str, help="GTF Annotation")
parser.add_argument('peaks_file', type=str, help="Peaks file from Homer")
parser.add_argument('--bed', action="store_const", const=True, default=False, 
	help="Input is BED format (default: homer)")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

