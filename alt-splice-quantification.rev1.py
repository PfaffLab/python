#!/usr/bin/env python
#
# 

# working out how to reduce a locus down to only unique alternative splicing features
# each feature will be an exon skipping event. so that means we need to know: 
# the exon that's skipped, 
# the junctions that include it, 
# the junction that skips it and 
# the exons adjacent to it that the skipping junction connects between

import sys,re,argparse
import subprocess as sp
from math import log
import numpy as np

def main(args):

	#
	# variables
	#

	if not check_for_required_stuff():
		sys.stderr.write("missing required tools...bailing out\n")
		return 1

	gtf = args.gtf
	reads = args.read_files
	genome_ref = args.genome_ref

	strand = ""
	temp_folder = "alt_splice.tmp"
	exons_temp = {}
	juncs_temp = {}
	exons_to_isoforms = {}
	junctions_to_exons = {}
	junctions_to_owner_exons = {}
	exons_to_owned_junctions = {}
	locus_rows = {}
	loci_to_names = {}
	# this is going to be a hash database of the alt locations. we'll use this to 
	# keep things organized and to be able to work with the pileup ouput after
	# aligning
	alocs = {}


	#
	# if the GTF file has the locus tag already we're good, otherwise we need to call on
	# gffread to add that tag in for us
	#
	if not gtf_has_locus_tag(gtf):
		sys.stderr.write("> locus tag was not found in your GTF. Using gffread to group transcripts into loci\n")
		fout = open(gtf + ".locus",'w')
		p1 = sp.Popen("gffread --cluster-only -F -T -o- {:s}".format(gtf).split(),stdout=fout)
		p1.wait()
		fout.close()
		gtf = gtf + ".locus"

	#
	# now we need to parse the GTF in so we have all of its information read to roll
	#
	sys.stderr.write("> parsing GTF and removing loci that contain only 1 transcript...\n")
	locus_rows,loci_to_names = gtf_to_hash(gtf)

	#
	# open connection to /dev/null to send unwanted output from system calls
	#
	fnull = open('/dev/null','w')
	al_fout = open("alt_locs.gtf","w")
	al_idx = 0

	#
	# process each locus to produce the alternative splice location GTF file
	#
	sys.stderr.write("> searching for alternative splice locations at each locus\n")
	for rloc in sorted(locus_rows.keys()):
		# initalize various hashes used to process each locus
		exons_temp = {}
		juncs_temp = {}
		exons_to_isoforms = {}
		junctions_to_exons = {}
		junctions_to_owner_exons = {}
		exons_to_owned_junctions = {}		

		#
		# make a temp GTF file so that we can process this locus using system commands
		#
		fout = open("temp_locus.gtf","w")
		for i in range(len(locus_rows[rloc])):
			fout.write(locus_rows[rloc][i] + "\n")
		fout.close()

		#
		# extract junctions using tophat's utility. this produces a bed format file with only
		# unique junctions. we'll pass these into a new file called temp.juncs
		#
		p1 = sp.Popen(["gtf_juncs","temp_locus.gtf"],stdout=sp.PIPE,stderr=fnull)
		# lines = p1.communicate()[0].split("\n")

		idx = 0
		# write junctions to file so we can use them with bedtools but also parse them
		# into a hash for quick access within the code
		fout = open("temp.juncs","w")

		for szl in p1.stdout:
			szl = szl.strip()
			if len(szl) > 0:
				arl = szl.split("\t")
				arl[2] = str(int(arl[2])+1)
				arl[3] = "JUNC{:04d}".format(idx)
				fout.write("\t".join(arl) + "\n")
				idx += 1
				juncs_temp[arl[3]] = arl[0:3]					

		fout.close()

		# strip out exons
		#fout = open(gtf + ".exons","w")
		#cmd = "grep exon " + gtf
		#p1 = sp.Popen(cmd.split(),stderr=fnull,stdout=fout)
		#fout.close()

		# convert to bed and collapse annotation to only unique exons
		fout = open("temp.exons","w")
		cmd = "gtf-to-bed temp_locus.gtf"
		p1 = sp.Popen(cmd.split(),stdout=sp.PIPE,stderr=fnull)
		p2 = sp.Popen("sort -k1,1 -k2,2n -k3,3n".split(),stdin=p1.stdout,stdout=sp.PIPE)
		p3 = sp.Popen("groupBy -i stdin -g 1,2,3 -c 4,6 -o distinct,distinct".split(),stdout=sp.PIPE,stdin=p2.stdout)
		#lines = p3.communicate()[0].split("\n")

		idx = 0
		for szl in p3.stdout:
			szl = szl.strip()
			if len(szl) > 0:
				arl = szl.split("\t")
				eid = "EXON{:04d}".format(idx)
				exons_to_isoforms[eid] = arl[3]
				arl[3] = eid
				fout.write("\t".join(arl) + "\n")
				idx += 1
				exons_temp[eid] = arl

		fout.close()
				
		# 
		# search for exons that appear to be skipped by junctions
		#

		# intersect the exons with the junctions checking for overlaps where the exon is 100% covered by the
		# junction. we'll make a hash of the junction id's containing the list of exons that the 
		# junction skips

		p1 = sp.Popen("bedtools intersect -wo -f 1 -a temp.exons -b temp.juncs".split(),stdout=sp.PIPE)
		#lines = p1.communicate()[0].split("\n")

		for szl in p1.stdout:
			szl = szl.strip()
			if len(szl) > 0:

				arl = szl.split("\t")

				# build a hash of the unique junction ids that are exon skipping junctions.
				# for each id we'll have a set of exon id's that are the exons skipped by 
				# the junction
				if arl[8] not in junctions_to_exons:
					junctions_to_exons[arl[8]] = set([])

				junctions_to_exons[arl[8]].update([arl[3]])

		# 
		# only continue if something was returned from the intersection
		#
		if len(junctions_to_exons.keys()) > 0:

			#
			# build two hashes that link junctions to exons and exons to junctions so we can look
			# up whatever we need.
			#

			#
			# NOTE: any single exon transcript, like a micro-rna, is going to show up as a skipped
			# exon but the exons_to_owned_junctions hash will contain no junctions for it.
			#

			p1 = sp.Popen("bedtools intersect -wo -a temp.exons -b temp.juncs".split(),stdout=sp.PIPE)
			#lines = p1.communicate()[0].split("\n")
			for szl in p1.stdout:
				szl = szl.strip()
				if len(szl) > 0:
					arl = szl.split("\t")
					
					jid = arl[8]
					eid = arl[3]

					if jid not in junctions_to_owner_exons:
						junctions_to_owner_exons[jid] = set([])

					if eid not in exons_to_owned_junctions:
						exons_to_owned_junctions[eid] = set([])

					if int(arl[9]) == 1:
						# check the exons already associated with the currend jid. we only want two exons
						# per junction otherwise things get messy. see comment following this section.
						
						# so what this loop will do is check the current exon against those already 
						# associated with the junction for overlaps. if they overlap then we'll just 
						# keep the one already in the set otherwise we'll add the new exon to the set
						found_overlap = False
						if len(junctions_to_owner_exons[jid]) > 0:
							ce_info = exons_temp[eid]
							for eid_j in list(junctions_to_owner_exons[jid]):
								e_info = exons_temp[eid_j]
								if ce_info[2] > e_info[1] and ce_info[1] < e_info[2]:
									found_overlap = True

						if not found_overlap:
							# ok to add this exon to the list
							junctions_to_owner_exons[jid].update([eid])

						exons_to_owned_junctions[eid].update([jid])

			#
			# we need to prune the information just collected so that each junction is only associated with two exons. sometimes
			# a junction will associate with multiple exons even though the locus has been reduced to a set of unique exons. this 
			# is because sometimes there are exons with the same start position but different lengths - like multiple 3' exons
			# that all have the same intron before them but have different lengths. I don't think it matters which exon we keep
			# just as long as the junction has two non-overlapping exons in its list
			#

			#
			# check through junctions_to_exons for any junction that says it is skipping an exon that
			# has no junctions (like a single exon transcript)
			#
			for eid in exons_to_owned_junctions.keys():
				if len(list(exons_to_owned_junctions[eid])) == 0:
					#
					# this exon has no junctions. make sure it is not involved in any of the 
					# downstream analysis
					#
					jids = junctions_to_exons.keys()
					for jid in jids:
						eset = list(junctions_to_exons[jid])
						elist = ",".join(eset)
						m = re.search(eid,elist)
						if m:
							# found it
							if len(eset) > 1:
								# junction hits other exons so we don't need to remove the entire junction
								# just this exon id
								tset = set(eset).difference(set([eid]))
								junctions_to_exons[jid] = list(tset)
							else:
								# there is only one exon associated with this junction so we can 
								# remove the entire event
								del junctions_to_exons[jid]

			#
			# loop through the junctions that skip exons and generate GTF output for each 
			# event
			#

			idx = 0
			tidx = 0
			for key in junctions_to_exons.keys():
				#
				# get exons that this junction connects - this is the skipping event
				#
				
				elist = list(junctions_to_owner_exons[key])				
				
				#
				# make sure these exons are sorted according to their genomic location
				#
				found_swap = True
				while found_swap:
					found_swap = False
					for i in range(len(elist)-1):
						if exons_temp[elist[i+1]][1] < exons_temp[elist[i]][1]:
							# swap
							temp = elist[i]
							elist[i] = elist[i+1]
							elist[i+1] = temp
							found_swap = True				

				alid = "ALOC_{:08d}".format(al_idx)

				# build information for the skip event
				alocs[alid] = {}
				alocs[alid]['loc'] = juncs_temp[key][0] + ":" + juncs_temp[key][1] + "-" + juncs_temp[key][2]
				alocs[alid]['locus'] = rloc
				alocs[alid]['skip'] = {}
				alocs[alid]['exons'] = {}
				alocs[alid]['skip']['elist'] = elist
				alocs[alid]['skip']['tset'] = set(exons_to_isoforms[elist[0]].split(",")).intersection(set(exons_to_isoforms[elist[1]].split(",")))
				alocs[alid]['skip']['splice_pos'] = pair_splice_position([exons_temp[elist[0]],exons_temp[elist[1]]])

				for eid in elist:
					e_info = exons_temp[eid]
					g_row = [
						e_info[0],
						"alt_splice",
						"exon",
						int(e_info[1])+1,
						e_info[2],
						"0.0",
						e_info[4],
						".",
						"gene_id \"{:s}\"; transcript_id \"{:s}\"; locus \"{:s}\";".format(alid,alid + ".skip",rloc)]

					al_fout.write("\t".join(map(str,g_row)) + "\n")

				#
				# now with the list of exons skipped by this junction we need to gather the junctions that
				# splice it in and then the exons that the junction belongs to
				#

				# get skipped exon(s)
				# tidx = 0
				eidx = 0
				elist = junctions_to_exons[key]
				for eid in elist:
					# get junctions that hook this exon in
					jlist = exons_to_owned_junctions[eid]
					# loop through these to make GTF rows for each pair of exons
					tidx = 0
					al_ei = "e" + str(eidx)

					alocs[alid]['exons'][al_ei] = []

					for jid in jlist:
						elist_sub = list(junctions_to_owner_exons[jid])
						
						# make sure the exons are sorted according to genomic location
						found_swap = True
						while found_swap:
							found_swap = False
							for i in range(len(elist_sub)-1):
								if exons_temp[elist_sub[i+1]][1] < exons_temp[elist_sub[i]][1]:
									# swap
									temp = elist_sub[i]
									elist_sub[i] = elist_sub[i+1]
									elist_sub[i+1] = temp
									found_swap = True				
						
						alocs[alid]['exons'][al_ei].append({})
						alocs[alid]['exons'][al_ei][tidx]['elist'] = elist_sub
						alocs[alid]['exons'][al_ei][tidx]['tset'] = set(exons_to_isoforms[elist_sub[0]].split(",")).intersection(set(exons_to_isoforms[elist_sub[1]].split(",")))
						alocs[alid]['exons'][al_ei][tidx]['splice_pos'] = pair_splice_position([exons_temp[elist_sub[0]],exons_temp[elist_sub[1]]])

						for eid_sub in elist_sub:
							e_info = exons_temp[eid_sub]
							g_row = [
								e_info[0],
								"alt_splice",
								"exon",
								int(e_info[1])+1,
								e_info[2],
								"0.0",
								e_info[4],
								".",
								"gene_id \"{:s}\"; transcript_id \"{:s}\"; locus \"{:s}\";".format(alid,alid + "." + al_ei + "." + str(tidx),rloc)]

							al_fout.write("\t".join(map(str,g_row)) + "\n")

						tidx += 1

					eidx += 1

				al_idx += 1

	al_fout.close()

	#
	# build fasta reference using gffread and the gtf created in the previous section
	#
	sys.stderr.write("> building fasta reference from alternative splice locations...\n")
	p1 = sp.Popen("gffread -g {:s} -w alt_locs.fa alt_locs.gtf".format(args.genome_ref).split(),stdout=fnull,stderr=fnull)
	p1.wait()
	#
	# build bowtie index from the fasta
	#
	sys.stderr.write("> building bowtie index...\n")
	p1 = sp.Popen("bowtie-build --offrate 1 alt_locs.fa alt_locs".split(),stdout=fnull,stderr=fnull)
	p1.wait()

	#
	# align reads - this section can be put into a loop to handle multiple read files
	#

	for rfile in args.read_files:

		rfile_parts = rfile.split("/")
		fout_stub = rfile_parts[len(rfile_parts)-1].split(".")[0]
		fout_name = fout_stub + ".alt_out"

		sys.stderr.write("> aligning {:s} with bowtie...\n".format(rfile))
		#bw_cmd = "bowtie -p {:d} -n 2 -e 40 -a -S alt_locs -".format(args.num_threads)
		#sys.stderr.write("cmd: " + bw_cmd + "\n")
		#p1 = sp.Popen(['gunzip', '-c', rfile],stdout=sp.PIPE)
		#p2 = sp.Popen(bw_cmd.split(),stdin=p1.stdout,stdout=sp.PIPE)
		#p3 = sp.Popen("samtools view -bS -F 0x4 -o temp.bam -".split(),stdin=p2.stdout)
		#p3.wait()
		res = sp.check_output("file {:s}".format(rfile),shell=True)
		m = re.search("gzip",res)
		if m:
			sp.call("bash -c 'bowtie --offrate 1 -p {:d} -a -n 1 -l 25 -e 60 -S alt_locs <(gunzip -c {:s}) | samtools view -bS -F 0x04 -o temp.bam -'".format(args.num_threads,rfile),shell=True)
		else:
			sp.call("bash -c 'bowtie --offrate 1 -p {:d} -a -n 1 -l 25 -e 60 -S alt_locs {:s} | samtools view -bS -F 0x04 -o temp.bam -'".format(args.num_threads,rfile),shell=True)

		sys.stderr.write("> sorting alignments...\n")
		p1 = sp.Popen("samtools sort temp.bam alt_locs".split())
		p1.wait()
		
		#
		# translate the alignments to genomic coordinates and pull out those that alignments that are
		# spliced.
		#
		
		sys.stderr.write("> parsing alignments and counting hits at splices\n")
		p1 = sp.Popen("samtools view alt_locs.bam".split(),stdout=sp.PIPE)
		
		#
		# parse
		#
		last_ref = ""
		current_ref = ""
		line_buffer = []

		aloc_counts = {}
		
		for szl in p1.stdout:
			szl = szl.strip()
			
			if len(szl) > 0:
				sam_row = szl.split("\t")
			
				current_ref = sam_row[2]
				
				if current_ref != last_ref:
					
					if len(line_buffer) > 0:
						#print "quantifying {:s}".format(last_ref)
	
						#
						# translate collected alignments to genomic positions
						#
						
						aloc_counts[last_ref] = 0
						
						# write to a temp file
						fout = open("temp.sam",'w')
						for i in range(len(line_buffer)):
							fout.write(line_buffer[i] + "\n")
						
						fout.close()
											
						# translate and read back in
						p2 = sp.Popen("tsam-to-gsam alt_locs.gtf temp.sam".split(),stdout=sp.PIPE)
						
						# parse
						for szll in p2.stdout:
							szll = szll.strip()
							if len(szll) > 0:
								arl = szll.split("\t")
								ctypes,clengths = parse_cigar(arl[5])
								if "N" in ctypes:
									# this is a junction alignment, check the anchors
									count_me = True
									for i in range(len(ctypes)):
										if ctypes[i] == "M" and clengths[i] < 4:
											# anchor is too short, dont' count this hit
											count_me = False
											#sys.stderr.write("short anchor alignment: dropping\n")
										#
									
									if count_me:
										# count this hit
										aloc_counts[last_ref] += 1
									#if
								#if
							#if
						#for
					
					#if
					
					# finished parsing last ref, now reset things so we can work on the next one
					line_buffer = []
				
				#if
				
			#if 
			
			# append current SAM line to the buffer
			last_ref = current_ref
			line_buffer.append(szl)
		
		#for
		
		# deal with last reference

		if len(line_buffer) > 0:
			#
			# translate collected alignments to genomic positions
			#
			
			aloc_counts[last_ref] = 0
			
			# write to a temp file
			fout = open("temp.sam",'w')
			for i in range(len(line_buffer)):
				fout.write(line_buffer[i] + "\n")
			
			fout.close()
								
			# translate and read back in
			p2 = sp.Popen("tsam-to-gsam alt_locs.gtf temp.sam".split(),stdout=sp.PIPE)
			
			# parse
			for szll in p2.stdout:
				szll = szll.strip()
				if len(szll) > 0:
					arl = szll.split("\t")
					ctypes,clengths = parse_cigar(arl[5])
					if "N" in ctypes:
						# this is a junction alignment, check the anchors
						count_me = True
						for i in range(len(ctypes)):
							if ctypes[i] == "M" and clengths[i] < 4:
								# anchor is too short, dont' count this hit
								count_me = False
							#
						
						if count_me:
							# count this hit
							aloc_counts[last_ref] += 1
						#if
					#if
				#if
			#for
		
		#if		
		
		#print aloc_counts
		
		if False:
		
			p1 = sp.Popen("samtools faidx alt_locs.fa".split())
			p1.wait()
	
			#
			# process the pileup data
			#
	
			sys.stderr.write("> processing pileup at splice positions\n")
			p1 = sp.Popen("samtools mpileup -f alt_locs.fa alt_locs.bam".split(),stdout=sp.PIPE,stderr=fnull)
	
			#fin = open("alt_locs.pileup","r")
			cid = ""
			ppos = 0
			aloc_counts = {}
			idl = []
	
			for szl in p1.stdout:
				szl = szl.strip()
				if len(szl) > 0:
					arl = szl.split("\t")
	
					if arl[0] != cid:
	
						if cid != "":
							aloc_counts[cid] = pcounts
	
						# 
						# new feature id
						idl = arl[0].split(".")
						
						# find the splice position for this alt loc id
						if idl[1] == "skip":
							# it's the skip event
							ppos = alocs[idl[0]]["skip"]["splice_pos"]
						else:
							# it's one of the skipped exon events
							ppos = alocs[idl[0]]['exons'][idl[1]][int(idl[2])]['splice_pos']
	
						# initalize counts
						pcounts = [0,0]
						# set id for the next round though the loop
						cid = arl[0]
	
					# if the position of the pileup data is at the position of the splice we
					# want to check the count of matches to the reference.
					if int(arl[1]) == ppos or int(arl[1]) == ppos+1:
						pdata = arl[4]
						# drop insertions
						pdata = re.sub("\+[0-9]+[ACGTNacgtn]+","",pdata)
						# drop deletions
						pdata = re.sub("\-[0-9]+[ACGTNacgtn]+","",pdata)
	
						i = 0
						while i < len(pdata):
							if pdata[i] == "$":
								# skip ends of reads
								i += 1
							elif pdata[i] == "^":
								# skip starts of reads
								i += 2
							elif pdata[i] == ',' or pdata[i] == '.':
								pcounts[ppos-int(arl[1])] += 1
	
							i += 1
	
			# insert the last one
			aloc_counts[cid] = pcounts
			#fin.close()
			
		#eif

		#
		# so now using the alocs strucure we can organize the counts information into a coherent output
		#

		sys.stderr.write("> reporting\n")
		al_out = open(fout_name,'w')

		# write header
		al_out.write("\t".join(["al_uid", "locus_id", "gene_name", "locus_tid", "aloc_id", "skip_tid", "skip_loc", "aloc_eid", "exon_tid", "skip_count", "inc_count", "skip_to_inc_ratio"]) + "\n")

		for key in sorted(alocs.keys()):
			#print alocs[key]
			# find skip event

			# look for the skip event for this alt location
			mean_count = 0
			al_key = key + ".skip"
			if al_key in aloc_counts:
				# skip_e = aloc_counts[al_key]
				# mean_count = (skip_e[0]+skip_e[1])/2.
				mean_count = aloc_counts[al_key]

			# the skip even could have skipped multiple exons so we need to generate output for 
			# each one and each one will have 1 or 2 supporting counts
			num_exons = len(alocs[key]['exons'])
			for i in range(num_exons):
				al_ei = "e" + str(i)
				# get number of junctions that support this exon
				num_juncs = len(alocs[key]['exons'][al_ei])
				p_counts = []
				# make a set of the transcript ids this exon is included in
				tlist = set([])
				for j in range(num_juncs):
					al_key = key + "." + al_ei + "." + str(j)
					if al_key in aloc_counts:
						p_counts += [aloc_counts[al_key]]
					else:
						p_counts += [0]

					if len(tlist) == 0:
						tlist = alocs[key]['exons'][al_ei][j]['tset']
					else:
						tlist.update(alocs[key]['exons'][al_ei][j]['tset'])

				e_mean_count = sum(p_counts)/(len(p_counts)*1.0)

				if False:
					the_ratio = 0
					if e_mean_count > 0 and mean_count == 0:
						the_ratio = -np.inf
					elif e_mean_count == 0 and mean_count > 0:
						the_ratio = np.inf
					elif e_mean_count > 0 and mean_count > 0:
						the_ratio = log2(mean_count)-log2(e_mean_count)
					else:
						the_ratio = 0

				the_ratio = log2(mean_count+1) - log2(e_mean_count+1)

				rloc = alocs[key]['locus']
				g_name = ",".join(list(loci_to_names[rloc][0]))
				t_ids = ",".join(list(loci_to_names[rloc][1]))

				# make a list of the transcript id(s) the skipping junction belongs to by subtracting out the
				# transcript id(s) of the current exon
				j_tset = alocs[key]['skip']['tset'].difference(tlist)

				al_out.write("\t".join(map(str,[
					key + "." + al_ei,
					rloc,
					g_name,
					t_ids,
					key,
					",".join(list(j_tset)),
					alocs[key]['loc'],
					al_ei,
					",".join(list(tlist)),
					mean_count,
					e_mean_count,
					the_ratio
				])) + "\n")

				#",".join(list(alocs[key]['skip']['tset'])),
				#",".join(list(tlist)),

		al_out.close()
		
		if args.keep_bam:
			sp.call("mv alt_locs.bam {:s}".format(fout_stub + ".bam"),shell=True)
		else:
			sp.call("rm alt_locs.bam",shell=True)
		
	#efor

	fnull.close()

	# clean up
	p1 = sp.Popen("rm temp.bam temp.sam temp.exons temp.juncs temp_locus.gtf".split())
	p1.wait()
	

	return 0

#def pileup_count(pdata):

def check_for_required_stuff():
	#
	# check for bowtie
	#
	try:
		res = sp.check_output("which bowtie",shell=True)
	except sp.CalledProcessError,e:
		print "bowtie is missing. Please install bowtie and make sure it is your system path"
		return False

	#
	# check for gffread
	#
	try:
		res = sp.check_output("which gffread",shell=True)
	except sp.CalledProcessError,e:
		print "gffread (via cufflinks) is missing. Please install cufflinks and make sure it is your system path"
		return False

	#
	# check for gtf_juncs
	#
	try:
		res = sp.check_output("which gtf_juncs",shell=True)
	except sp.CalledProcessError,e:
		print "gtf_juncs (via tophat) is missing. Please install tophat and make sure it is your system path"
		return False

	#
	# check for gtf-to-bed
	#
	try:
		res = sp.check_output("which gtf-to-bed",shell=True)
	except sp.CalledProcessError,e:
		print "gtf-to-bed (via Driscoll) is missing. Make sure it's accessible via the system path."
		return False

	#
	# check for tsam-to-gsam
	#
	try:
		res = sp.check_output("which tsam-to-gsam",shell=True)
	except sp.CalledProcessError,e:
		print "tsam-to-gsam (via Driscoll) is missing. Make sure it's accessible via the system path."
		return False


	return True


def gtf_has_locus_tag(gtf):

	#
	# check if the GTF file contains the locus tag
	#
	fin = open(gtf,'r')
	szl = fin.readline()
	fin.close()

	m = re.search('locus',szl)
	if m:
		return True

	return False


def parse_cigar(sz):

	split_1 = re.split("[A-Z]",sz)
	split_2 = re.split("[0-9]+",sz)

	split_1 = split_1[0:(len(split_1)-1)]
	split_2 = split_2[1:]

	return (split_2,map(int,split_1))

def parse_gtf_attr(field):
	#
	# parse the attributes field of a gtf row into a hash
	#
	fsplit = field.split("\"")
	attrs = {}

	n = len(fsplit)-1
	i = 0
	while i < n:
		key = re.sub(';','',fsplit[i])
		attrs[key.strip()] = fsplit[i+1].strip()
		i += 2

	return attrs

def gtf_to_bed(row):
	#
	# reformat a gtf row into a bed row
	#
	fsplit = row.split("\t")
	attrs = parse_gtf_attr(fsplit[8])

	lout = [
		fsplit[0],
		str(int(fsplit[3])-1),
		fsplit[4],
		attrs['transcript_id'],
		'0.0',
		fsplit[6]
	]

	return "\t".join(lout)

def gtf_to_hash(gtf):
	#
	# parse the GTF information into a hash by locus. each hash position will hold a list
	# of the raw rows from the GTF. we also want a table of locus to gene name and transcript
	# id for the final output
	#

	parsed = {}
	name_table = {}
	fin = open(gtf,'r')

	for szl in fin:
		szl = szl.strip()
		if len(szl) > 0:
			arl = szl.split("\t")
			if arl[2] == 'exon':
				attrs = parse_gtf_attr(arl[8])
				if attrs['locus'] not in parsed:
					parsed[attrs['locus']] = []
					name_table[attrs['locus']] = [set([]),set([])]
				
				parsed[attrs['locus']].append(szl)
				name_table[attrs['locus']][0].update([attrs['gene_id']])
				name_table[attrs['locus']][1].update([attrs['transcript_id']])

	fin.close()

	rlocs = name_table.keys()

	# remove loci that contain single transcripts
	for key in rlocs:
		if len(list(name_table[key][1])) == 1:
			del name_table[key]
			del parsed[key]

	return (parsed,name_table)



def pair_splice_position(elist):
	#
	# return the position within a continuous exon pair sequence of the junction
	# between the two exons.
	#

	# check strand
	if elist[0][4] == '-':
		# get length of second feature
		return int(elist[1][2])-int(elist[1][1])
	
	# get length of the first feature
	return int(elist[0][2])-int(elist[0][1])

def log2(x):
	if x == 0:
		return np.inf

	return log(x)/log(2)

#
# entry point
#
parser = argparse.ArgumentParser(description="Alternative splicing pipeline.")
parser.add_argument('gtf',type=str,help="GTF annotation")
parser.add_argument('genome_ref',type=str,help="FASTA Genome reference (index required) from which sequence based on the GTF annotation can be extracted")
parser.add_argument('read_files',type=str,metavar='reads',nargs='+',help='FASTQ reads to process against the alternative splice database')
parser.add_argument('--old-quals',dest='old_quals',default=False,action='store_const',const=True,help='Reads use phred+64 (solexa 1.3) qualities (default: False)')
parser.add_argument('--keep-bam', dest='keep_bam', default=False, action='store_const', const=True, help="Keep the BAM alignments to the alternative splice locations (default: off)")
parser.add_argument('-p',dest="num_threads",default=1,type=int,help="Number of processors to use for bowtie stage (default: 1)")
args = parser.parse_args()

if __name__ == "__main__":
	sys.exit(main(args))
