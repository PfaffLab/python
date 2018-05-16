#!/usr/bin/python
#
# Shawn Driscoll
# 20160924
#

#==============================================================================
# imports
#==============================================================================

import sys, argparse, os, re
import hashlib
import gzip

#==============================================================================
# globals
#==============================================================================

#==============================================================================
# main
#==============================================================================

def main(argv):

	# variables
	rres = None

	#
	# we need a list of all of the fastq.gz files within the current folder
	rres = os.popen("find `pwd` -name \"*.fastq.gz\"").readlines()
	rres += os.popen("find `pwd` -name \"*.fq.gz\"").readlines()

	#
	# loop through rres and sort the samples into a dict. 
	dfirst = {}
	for i in range(len(rres)):
		#
		# clean and split up the full path to the fastq.gz file
		tmp = rres[i].strip().split("/")
		#
		# sample id is the folder name containing the file
		sid = tmp[-2]
		#
		# file name is the last part
		fname = tmp[-1]
		#
		# if the sample isn't in the dict yet then add a slot for it
		if sid not in dfirst:
			ll = len(tmp)
			bpath = "/".join(tmp[0:(ll-1)]) + "/"
			dfirst[sid] = dict(left=[], right=[], basepath=bpath)

		#
		# check if the current file is a left or right side read
		r = re.search("[^A-Z0-9]R2[^A-Z0-9]", fname)
		if r:
			# right side
			dfirst[sid]['right'].append(fname)
		else:
			# left side
			dfirst[sid]['left'].append(fname)


	#
	# done sorting. print out a tab-delim table of the samples, the base path
	# and the left and right reads.
	#

	#
	# header
	print "sample\tbasepath\tleft_reads\tright_reads"
	#
	# loop through
	for sid in sorted(dfirst.keys()):
		lout = [sid, dfirst[sid]['basepath'], ",".join(sorted(dfirst[sid]['left'])), 
			",".join(sorted(dfirst[sid]['right']))]
		print "\t".join(lout)



	return 0


##
# in this one we  use the file names to establish samples
def main2(argv):

	# variables
	rres = None

	#
	# we need a list of all of the fastq.gz files within the current folder
	rres = os.popen("find `pwd` -name \"*.fastq.gz\"").readlines()
	rres += os.popen("find `pwd` -name \"*.fq.gz\"").readlines()

	#
	# loop through rres and sort the samples into a dict. 
	dfirst = {}


	lsamples = []
	hhit = set()

	for p in rres:
		p = p.strip()
		pparts = bustup(p)
		if pparts[3]:
			# first mate file
			hh = hashlib.md5(p).hexdigest()
			dfirst[hh] = pparts

	# loop back through to pair stuff up
	for p in rres:
		p = p.strip()
		pparts = bustup(p)

		if not pparts[3]:
			# second mate file
			phat = re.sub("\_R2\_", "_R1_", p)
			hh = hashlib.md5(phat).hexdigest()
			if hh in dfirst:
				# found a pair
				bpath = pparts[0]
				sid = re.sub("\_R2\_.+", "", pparts[1])
				mate1 = dfirst[hh][1]
				mate2 = pparts[1]
				rtype = "PE-{}".format(pparts[4])
				if argv.just_reads:
					lsamples.append(["{}/{}".format(bpath, mate1), "{}/{}".format(bpath, mate2)])
				else:
				# make a hollow version of the output I get in the samples excel files
				# type,geno,age,sort,labber,location,year,read-type,stranded,quals,sample_ref,source_path,files
					lout = ["u", "u", "u", "u", "AUTO", "MARS", "1999Dec", rtype, "R", "33", re.sub("_", ".", sid), bpath, mate1]
					if mate2 != "none":
						lout.append(mate2)
					#lsamples.append([sid, bpath, rtype, mate1, mate2])
					lsamples.append(list(lout))
					#lsamples.append([sid, bpath, rtype, mate1, mate2])
				hhit.add(hh)

	#
	# confirm that we got them all
	for hh in dfirst.keys():
		if hh not in hhit:
			# we have a single-end case
			bpath = dfirst[hh][0]
			sid = re.sub("\_R1\_.+", "", dfirst[hh][1])
			mate1 = dfirst[hh][1]
			mate2 = "none"
			rtype = "SE-{}".format(dfirst[hh][4])
			if argv.just_reads:
				lsamples.append(["{}/{}".format(bpath, mate1), "{}/{}".format(bpath, mate2)])
			else:
				# make a hollow version of the output I get in the samples excel files
				# type,geno,age,sort,labber,location,year,read-type,stranded,quals,sample_ref,source_path,files
				lout = ["u", "u", "u", "u", "AUTO", "MARS", "1999Dec", rtype, "R", "33", sid, bpath, mate1]
				if mate2 != "none":
					lout.append(mate2)
				#lsamples.append([sid, bpath, rtype, mate1, mate2])
				lsamples.append(list(lout))

	if not args.just_reads:
		print "type\tgeno\tage\tsort\tlabber\tlocation\tyear\tread_type\tstranded\tquals\tsample_ref\tsource_path\tfiles\tfiles2"
		
	for l in lsamples:
		print "\t".join(l)

	return 0


# provide full path of file this will break it apart a little
def bustup(p):

	pparts = p.strip().split("/")
	n = len(pparts)
	fname = pparts[n-1]
	
	# if path started with a '/' then the broken version will have a blank as the 
	# first list element.
	if p[0] == "/":
		bpath = "/" + "/".join(pparts[1:(n-1)])
	else:
		bpath = "/".join(pparts[0:(n-1)])

	fparts = fname.split(".")
	fstub = fparts[0]

	first_mate = True
	if re.search("\_R2\_", fstub):
		first_mate = False

	# get read length
	rlen = get_read_length(p)

	lout = [bpath, fname, fstub, first_mate, rlen]
	return lout


def get_read_length(f):

	rlen = 0

	##
	## if file is gzipped then we have to call the system to deal with it
	if re.search("\.gz$", f):

		with gzip.open(f, "r") as fin:
			szl = fin.readline()
			szl = fin.readline().strip()
			rlen = len(szl)

	else:

		with open(f, "r") as fin:
			szl = fin.readline()
			szl = fin.readline().strip()
			rlen = len(szl)

	return rlen


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="About.")
#parser.add_argument('fin', type=str, help="Input file")

parser.add_argument('-l', type=int, default=0, action="store", 
	help="This value indicates the number of folder levels to keep connected to the file name")

parser.add_argument('--just-reads', action="store_const", const=True, default=False, 
	help="Just output a two-column list of the reads with the full path attached.")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main2(args))
		#sys.exit(main())
	except KeyboardInterrupt:
		sys.stderr.write("\n\nUser killed execution...\n\n")
