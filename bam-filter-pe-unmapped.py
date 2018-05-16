#!/usr/bin/env python
#
# Shawn Driscoll
# 20130320
#
# Filter out unmapped reads from paired end alignments. Alignments should be
# sorted by name (samtools sort -n).
#

import sys, argparse
import subprocess as sp

#==============================================================================
# main
#==============================================================================
def main(args):

	# confirm data is paired


	# get the sam header
	p1 = sp.Popen("samtools view -H {:s}".format(args.bam).split(), stdout=sp.PIPE)
	sam_header = ""
	for szl in p1.stdout:
		sam_header += szl

#	p1.stdout.close()

	# open output stream

#	if args.ofile == "stdout":
#		p1 = sp.Popen("samtools view -S -".split(), stdin=sp.PIPE, stdout=sys.stdout)
#	else:
#		p1 = sp.Popen("samtools view -S -o {:s} -".format(args.ofile).split(), stdin=sp.PIPE)

	# open stream from bam
	p2 = sp.Popen("samtools view {:s}".format(args.bam).split(), stdout=sp.PIPE)

	# write header
	sys.stdout.write(sam_header)
#	p1.stdin.write(sam_header)
	
	lbuffer = []

	for szl in p2.stdout:
		ll = szl.strip().split("\t")
		
		# ignore secondary alignments
		if not int(ll[1]) & 0x100:
			# append current line to buffer
			lbuffer.append(list(ll))

		if len(lbuffer) == 2:
			flag1 = int(lbuffer[0][1])
			flag2 = int(lbuffer[1][1])

			# compare names - might need to trim off the /1 or /2
			if lbuffer[0][0] == lbuffer[1][0]:
				# names match
				if (flag1 & 0x40 and flag2 & 0x80) or (flag1 & 0x80 and flag2 & 0x40):
					# these are a pair
					if flag1 & 0x4 and flag2 & 0x4:
						# one or both are unaligned. adjust flag so both are set to unaligned
						if not flag1 & 0x4:
							flag1 += 4
						if not flag2 & 0x4:
							flag2 += 4

						lbuffer[0][1] = str(flag1)
						lbuffer[1][1] = str(flag2)

						# print back out
#						p1.stdin.write("\t".join(lbuffer[0]) + "\n")
#						p1.stdin.write("\t".join(lbuffer[1]) + "\n")
						sys.stdout.write("\t".join(lbuffer[0]) + "\n")
						sys.stdout.write("\t".join(lbuffer[1]) + "\n")

				else:
					sys.stderr.write("Error: neighboring alignments are not part of pair!\n")
					sys.stderr.write("\t".join(lbuffer[0]) + "\n")
					sys.stderr.write("\t".join(lbuffer[1]) + "\n")
#					p1.stdin.close()
					return(1)

			else:
				sys.stderr.write("Error: neighboring alignments are not part of pair!\n")
#				p1.stdin.close()
				return(1)

			lbuffer = []


#	p1.stdin.close()


#==============================================================================
# main entry point
#==============================================================================

parser = argparse.ArgumentParser(description="Filter out unmapped reads from paired end alignments. Alignments should be sorted by name (samtools sort -n).")
parser.add_argument("bam", type=str, help="BAM file")
parser.add_argument("-o", dest="ofile", type=str, default="stdout", help="Output file (default: stdout)")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
