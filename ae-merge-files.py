#!/usr/bin/env python
#
# ae-merge-files.py
#
# Shawn Driscoll
# 20130128
#
# Merges the output of 2 or more ae quantification files (from ae-quantify-from-bam)
#

import sys, argparse

#
# MAIN
#
def main(args):

	num_files = len(args.ae_files)

	if num_files < 2:
		sys.stderr.write("No need to merge one file with itself, broseph.\n")
		return(1)

	ae_hash = {}
	ae_id = ""
	col_stub = ""

	
	# parse first file

	fin = open(args.ae_files[0], "r")
	# snag the header
	ll = fin.readline().strip().split("\t")
	header = list(ll[0:9])
	# set count column names
	col_stub = args.ae_files[0].split(".")[0]
	
#	header.append(col_stub + "_inc1")
#	header.append(col_stub + "_inc2")
#	header.append(col_stub + "_skip")

	header.append(col_stub + "_depth")
	header.append(col_stub + "_pic")

	for szl in fin:
		ll = szl.strip().split("\t")
		inc = max(float(ll[9]), float(ll[10]))
		skip = float(ll[11])
		depth = inc+skip

		if depth > 0:
			pic = inc/depth
		else:
			pic = 0

		ll_new = list(ll[0:9]) + [depth, pic]

		ae_hash[ll[0]] = list(ll_new)

	fin.close()

	# now loop through the rest of the files and append their counts
	for i in range(1,len(args.ae_files)):
		# deal with header additions for this file
		col_stub = args.ae_files[i].split(".")[0]
#		header.append(col_stub + "_inc1")
#		header.append(col_stub + "_inc2")
#		header.append(col_stub + "_skip")
		header.append(col_stub + "_depth")
		header.append(col_stub + "_pic")

		# open file and parse in its counts
		fin = open(args.ae_files[i], "r")
		szl = fin.readline()
		for szl in fin:
			ll = szl.strip().split("\t")

			inc = max(float(ll[9]), float(ll[10]))
			skip = float(ll[11])
			depth = inc+skip
			if depth > 0:
				pic = inc/depth
			else:
				pic = 0

			# append counts to ae loc row
			ae_hash[ll[0]] += [depth, pic]

		fin.close()

	# now that we have it all parsed we can print it all back out

	# print header
	print "\t".join(header)

	for ae_id in sorted(ae_hash.keys()):
		print "\t".join(map(str, ae_hash[ae_id]))

#
# main entry point
#

parser = argparse.ArgumentParser(description="Merge the output from 2 or more AE quantification files.")
parser.add_argument('ae_files',type=str,metavar='ae',nargs='+',help='AE quantification file(s)')
# options

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
	