#!/usr/bin/python
#
# juncs-and-boundaries.py
# Shawn Driscoll
# 20140619
#
# Extract junctions and boundary counts (reads that cross junctions but 
# are not spliced) from sorted BAM alignments. 
#

import sys, argparse, pysam

# -----------------------------------------------------------------------------
# main and other functions
# -----------------------------------------------------------------------------

def main(args):
	# variables
	# buffers for each "locus" of alignments
	juncs = {}
	da = {}
	alns = []
	aln_last = [-1, -1, -1]
	gap = 0
	locus_idx = 0
	min_anchor = 0
	min_overlap = 0
	rlen = 0

	# take a quick peak to find read length and to check if the alignments
	# are coordinate sorted.
	if args.S:
		sin = pysam.Samfile(args.alignments, "r")
	else:
		sin = pysam.Samfile(args.alignments, "rb")

	sorted = True
	lasta = None
	nread = 0
	for aln in sin:
		if rlen == 0:
			rlen = len(aln.seq)
		
		if aln.flag & 0x4:
			continue
		
		nread += 1
		
		if lasta is not None:
			if aln.tid != lasta.tid:
				# different references. if the difference in the reference
				# index is less than zero then the reads aren't sorted
				if aln.tid-lasta.tid < 0:
					sorted = False
					break
			else:
				# same reference. if the current alignment's position is less
				# that the last one then these are not sorted.
				if aln.pos - lasta.pos < 0:
					sorted = False
					break
		
		lasta = aln
		if nread > 500:
			break
	
	sin.close()

	if not sorted:
		# kill it
		sys.stderr.write("Alignments do not appear to be coordinate sorted!\n")
		return(1)

	# compute anchor and overlap mins from the read length
	min_anchor = round(args.a * rlen)
	min_overlap = round(args.o * rlen)

	# open file
	if args.S:
		sin = pysam.Samfile(args.alignments, "r")
	else:
		sin = pysam.Samfile(args.alignments, "rb")
	
	# print header
	print "\t".join(["chrom", "start", "end", "count", "locus", "donor_ovl", "donor_use", "acc_ovl", "acc_use", "psi5", "psi3", "theta5", "theta3"])
	
	# read on
	for aln in sin:
		if aln.flag & 0x4 or aln.mapq < args.q:
			continue
		
		# first determine if we are in the same locus
		# if aln doesn't overlap aln_last and the gap is > args.g then
		# process the buffered reads verse the buffered donors and 
		# acceptors
		
		if aln_last[0] >= 0:
			gap = check_aln_overlap(aln_last, aln)
			if gap < 0 or gap > args.g:
				# process buffer
				if len(da.keys()) > 0:
					result = process_buffer(alns, da, min_anchor, min_overlap)
					locus_idx += 1
					# pass these into the da hash
					for k in result.keys():
						result[k][4] = locus_idx
						da[k] = list(result[k])
					
					print_results(juncs, da, sin.header)
					
				juncs = {}
				alns = []
				da = {}
				aln_last = [-1, -1, -1]
		
		# determine if read is spliced and find it's aligned extent
		spliced = False
		aleft = aln.pos+1
		aright = aleft
		for c in aln.cigar:
			if c[0] == 0 or c[0] == 2:
				aright += c[1]
			if c[0] == 3:
				spliced = True
				aright += c[1]
				
		
		if spliced:
			# get coordinates of junction and find overall left and right 
			# coordinates
			left = aln.pos+1
			left0 = left
			
			n_cigar = len(aln.cigar)
			
			for c in aln.cigar:
				if c[0] == 3:
					# left has been adjusted up to this position so make right side
					# and we have it 
					right = left+c[1]-1
					
					# make junction name
					jid = "{}:{}:{}".format(aln.tid, left, right)
					
					# check if it exists already. if not add it.
					if jid not in juncs:
						juncs[jid] = [aln.tid, left, right, 0]
					
					# increment count
					juncs[jid][3] += 1
					
					# get donor and acceptor ids and add them to a hash if they don't
					# exist
					donor_id = "{}:{}".format(aln.tid, left)
					if donor_id not in da:
						# these are: tid, left, donor/acceptor, overlap count, locusid, use count 
						da[donor_id] = [aln.tid, left, 0, 0, 0, 0]
						
					# increment usage count
					da[donor_id][5] += 1
					
					acc_id = "{}:{}".format(aln.tid, right)
					if acc_id not in da:
						da[acc_id] = [aln.tid, right, 1, 0, 0, 0]
					
					# increment usage count
					da[acc_id][5] += 1
					
					# continue on to find more splices if any are there
					left += c[1]
				elif c[0] == 0 or c[0] == 2:
					left += c[1]
				
			this_aln = [left0, left]
			
		else:
			# not a spliced alignment. find left and right coordinates and buffer it
#			left = aln.pos+1
#			right = left
#			for c in aln.cigar:
#				if c[0] == 0 or c[0] == 2 or c[0] == 3:
#					right += c[1]			
			alns.append([aleft, aright])
			this_aln = [aleft, aright]
		
		# set last boundary. don't update range if this one doesn't 
		# extend at all to the right.
		if not (aln_last[0] >= 0):
			aln_last = [aln.tid, this_aln[0], this_aln[1]]
		else:
			if this_aln[1] > aln_last[2]:
				aln_last[2] = this_aln[1]
	
	
	# handle the final bit
	if len(da.keys()) > 0:
		result = process_buffer(alns, da, min_anchor, min_overlap)
		locus_idx += 1
		# pass these into the da hash
		for k in result.keys():
			result[k][4] = locus_idx
			da[k] = list(result[k])
		
	# print this business out
	print_results(juncs, da, sin.header)

	sin.close()
	
	return(0)

# --
# print_results
# print out the junctions line by line with the donor/acceptor
# overlap hits and junction counts. handles a single locus at a time.
def print_results(j, d, header):
	# variables
	jout= []
	jlist = []
	
	# parse the junctions out from the hash and build a list
	for k in j.keys():
		jout = list(j[k])
		# get donor and acceptor
		id = "{}:{}".format(jout[0], jout[1])
		if id in d:
			donor = list(d[id])
		
		id = "{}:{}".format(jout[0], jout[2])
		if id in d:
			acc = list(d[id])
		
		# append locus index
		jout.append(donor[4])
		
		# append counts
		jout.append(donor[3]) # overlap
		jout.append(donor[5]) # toal usage
		jout.append(acc[3])   # overlap
		jout.append(acc[5])   # total usage
		
		# append splicing metrics
		jout.append(jout[3]*1.0/donor[5]) # psi5
		jout.append(jout[3]*1.0/acc[5])   # psi3
		jout.append(donor[5]*1.0/(donor[5]+donor[3])) # theta5
		jout.append(acc[5]*1.0/(acc[5]+acc[3]))       # theta3
		
		#print "\t".join(map(str, jout))
		jlist.append(list(jout))
	
	# sort junctions by left position
	jlist.sort(key=lambda x: x[1])
	
	# print the junctions out. as we go replace the 
	# rnames with the real ones in the header of the Samfile
	# object
	for jout in jlist:
		tmp = header['SQ'][jout[0]]['SN']
		jout[0] = tmp
		tmp = "JLOC_{:08d}".format(jout[4])
		jout[4] = tmp
		print "\t".join(map(str, jout))
		

# --
# check_aln_overlap
# Returns 0 if the two alignments overlap otherwise
# it is assumed that a2 is downstream of a1 and the 
# length of the gap between the is returned. If a1 and
# a2 are not aligned to the same reference then -1 is 
# returned.
# a1 is just a three value list with tid and left and right coordinates.
def check_aln_overlap(a1, a2):
	
	if a1[0] != a2.tid:
		return -1
	
	# get left and right of a1 (excluding soft clips)
	a1_left = a1[1]
	a1_right = a1[2]
	
	# get left and right of a2 (excluding soft clips)
	a2_left = a2.pos + 1
	a2_right = a2_left
	for c in a2.cigar:
		if c[0] == 0 or c[0] == 2 or c[0] == 3:
			a2_right += c[1]
	
	if a1_left < a2_right and a1_right > a2_left:
		# ovelap!
		return 0
	
	# no overlap, return lenght of gap
	return a2_left - a1_right + 1

# --
# process_buffer
# checks buffered alignments against the buffered donor and
# acceptor sites to count overlaps.
def process_buffer(alns, das, min_anchor, min_overlap):
	
	# loop through alns
	for r in alns:
		for k in das.keys():
			da = das[k]
			if da[2] == 0:
				# donor site.
				if (da[1] - r[0]) >= min_anchor and (r[1] - da[1]) >= min_overlap:
					# got it
					das[k][3] += 1
					
			else:
				# acceptor site.
				if (da[1] - r[0]) >= min_overlap and (r[1] - da[1]) >= min_anchor:
					# got it
					das[k][3] += 1
	
	return das

# -----------------------------------------------------------------------------
# entry point
# -----------------------------------------------------------------------------


parser = argparse.ArgumentParser(description="Extracts junctions and donor/acceptor overlaps from sorted BAM alignments. psi and theta splicing metrics are also calculated.")
parser.add_argument('alignments', type=str, help="BAM alignments (or SAM with -S)")
parser.add_argument("-S", dest="S", action="store_const", const=True, default=False, 
				help="Alignments are SAM (default expected: BAM)")
parser.add_argument("-q", dest="q", type=int, action="store", default=1, 
				help="Minimum MAPQ (default: 1)")
parser.add_argument("-g", dest="g", type=int, action="store", default=20,
				help="Maximum tolerated coverage gap before starting new locus (default: 20)")
parser.add_argument("-o", dest="o", type=float, action="store", default=0.1,
				help="Minimum overlap into intron for overlap hits as a ratio of read length (default: 0.1)")
parser.add_argument("-a", dest="a", type=float, action="store", default=0.04,
				help="Minimum anchor for junctions (default: 0.04)")

args = parser.parse_args()

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("killed it\n")
