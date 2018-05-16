#!/usr/bin/env python
#
# ae-quantify-from-bam.py
#
# Shawn Driscoll
# 20130127
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Quantify AE from a BAM file. Looks for the <name>.gtf.ae.bed file that would have
# been generated from make-alt-exon-bed.
#

import sys,argparse
import subprocess as sp

#
# main
#

def main(args):

	#
	# variables
	#
	out_file = args.bam + ".aehits"
	if len(args.stub) > 0:
		out_file = args.stub + ".aehits"

	#
	# check for commands
	#
	sys.stderr.write("checking for bedtools...\n")
	res = sp.call("which bedtools", shell=True, stderr=sp.STDOUT)
	if res == 1:
		sys.stderr.write("Error: bedtools is not in the system path\n")
		return 1
	
	# check for the AE bed file
	try:
		fin = open(args.gtf + ".ae.bed")
	except:
		sys.stderr.write("Couldn't find BED file for your AE file. I'll make one now...\n")
		cmd = "make-alt-exon-bed {:s}".format(args.gtf)
		p1 = sp.Popen(cmd.split())
		p1.wait()
	else:
		fin.close()
	
	# we need to parse the AE file into a hash
	ae_file = parse_ae_file(args.gtf + ".ae")
	ae_count_table = {}
	
	# build count table
	for aeid in sorted(ae_file.keys()):
		# counts will be inc1, inc2, skip
		ae_count_table[aeid] = [0, 0, 0]
	
	# now we can intersect the aligned reads with the bed file. we'll read the result of the 
	# intersection in a loop
	cmd1 = "bedtools bamtobed -split -i {:s}".format(args.bam)
	cmd2 = "bedtools intersect -wo -a stdin -b {:s}".format(args.gtf + ".ae.bed")
	
	# 10 = hits
	# 9 = ae info
	
	p1 = sp.Popen(cmd1.split(), stdout=sp.PIPE)
	p2 = sp.Popen(cmd2.split(), stdin=p1.stdout, stdout=sp.PIPE)
	
	for szl in p2.stdout:
		ll = szl.strip().split("\t")
		
		# first make sure that the hit length is equal to the ae location length
		
		ae_len = int(ll[8]) - int(ll[7])
		if ae_len == int(ll[10]):
			# good to go, entire region was covered by this alignment
			
			# split the name field
			ae_list = ll[9].split(",")
			for i in range(len(ae_list)):
				# AELOC_00000007.skip|uc011wht.1|0|intron	
				ae_loc = ae_list[i].split("|")
				
				# split the first field to get the AE location id and the sub id (inc1, inc2, skip)
				temp = ae_loc[0].split(".")
				ae_id = temp[0]
				ae_type = temp[1]
				
				# add the hit to the appropriate field
				if ae_type == "skip":
					ae_count_table[ae_id][2] += 1
				elif ae_type == "inc1":
					ae_count_table[ae_id][0] += 1
				elif ae_type == "inc2":
					ae_count_table[ae_id][1] += 1

	# format of AE file
#	ae_id	ae_genoName	ae_genoStart	ae_genoEnd	ae_transcript_ids	ae_gene_names	ae_exon_index	strand	exon_type	9:ae_left	10:ae_right	11:skip_genoName	12:skip_genoStart	13:skip_genoEnd	skip_transcript_ids	skip_gene_names	skip_intron_index	skip_left	skip_right

	
	# done counting hits, now we can build the final output table and write to file
	fout = open(out_file, "w")
	
	# write the header
	fout.write("ae_id\tgene_name\texon_transcript_id\texon_index\texon_type\texon_location\tskip_transcript_id\tintron_index\tintron_location\tinc1\tinc2\tskip\tpic\n")
	
	for ae_id in sorted(ae_file.keys()):
		ll = ae_file[ae_id]
		lcounts = ae_count_table[ae_id]
		l_out = []
		# ae id
		l_out.append(ae_id)
		
		#
		# append ae exon info
		#
		
		# gene
		l_out.append(ll[5])
		# transcript
		l_out.append(ll[4])
		# exon index
		l_out.append(ll[6])
		# exon type
		l_out.append(ll[8])
		# location
		l_out.append(ll[1] + ":" + ll[2] + "-" + ll[3])

		#
		# append skip info
		#

		# transcript
		l_out.append(ll[14])
		# intron index
		l_out.append(ll[16])
		# location
		l_out.append(ll[11] + ":" + ll[12] + "-" + ll[13])
		
		# 
		# append hits
		#
		l_out += list(lcounts)
		
		inc_tmp = max(lcounts[0], lcounts[1])
		skip_tmp = lcounts[2]
		pic_tmp = 0
		if inc_tmp+skip_tmp > 0:
			pic_tmp = inc_tmp/float(inc_tmp+skip_tmp)
		
		l_out.append("%0.4f" % pic_tmp)
		
		fout.write("\t".join(map(str, l_out)) + "\n")
	
	# done
	fout.close()


def parse_ae_file(fname):
		
	ae_table = {}	
	fin = open(fname, "r")
	
	# skip header
	szl = fin.readline()
	
	# loop through file and parse it in
	for szl in fin:
		ll = szl.strip().split("\t")
		ae_table[ll[0]] = list(ll)
	
	fin.close()
	return(ae_table)

#==============================================================================
# entry point
#==============================================================================

# modes list

# bowtie-express
# bowtie2-express
# bwa-express
# tophat-cufflinks

parser = argparse.ArgumentParser(description="g3nom1c!")
parser.add_argument("gtf", type=str, action="store", help="GTF file that matches the transcriptome your reads are aligned to")
parser.add_argument("bam", type=str, action="store", help="BAM alignments to a transcriptome.")
parser.add_argument("-o", type=str, dest="stub", action="store", default="", help="Output stub for files. Default is to use the BAM file name.")

args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
