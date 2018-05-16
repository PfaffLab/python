#!/usr/bin/python
#
# gtf-to-donor-acceptor-annotation.py
# Shawn Driscoll
# 20140619
#
# Extract junctions and boundary counts (reads that cross junctions but 
# are not spliced) from sorted BAM alignments. 
#

import sys, argparse

# -----------------------------------------------------------------------------
# main and other functions
# -----------------------------------------------------------------------------

def main(args):
	# variables
	dgtf = {}
	donors = {}
	acceptors = {}
	junctions = {}

	# load the gtf
	sys.stderr.write("parsing gtf into transcripts\n")
	dgtf = load_gtf(args.gtf)
	
	sys.stderr.write("creating junctions from exon boundaries\n")
	
	for k in dgtf.keys():
		if len(dgtf[k]["exons"]) == 1:
			# single exon transcripts don't have junctions!
			continue
		
		# make junctions
		n = len(dgtf[k]["exons"])
		for i in range(n):
			
			if i < (n-1):
				# end point is a donor (ignoring strand)
				did = "{:s}:{:d}".format(dgtf[k]["exons"][i][0], dgtf[k]["exons"][i][2]+1)
				
				if did not in donors:
					donors[did] = {"gene_id":[], "transcript_id":[], "gene_name":[], "strand":[]}
				
				donors[did]["gene_id"].append(dgtf[k]["info"]["gene_id"])
				donors[did]["transcript_id"].append(dgtf[k]["info"]["transcript_id"])
				donors[did]["gene_name"].append(dgtf[k]["info"]["gene_name"])
				donors[did]["strand"].append(dgtf[k]["info"]["strand"])
			
			if i > 0:
				# start is an acceptor (ignoring strand)
				aid = "{:s}:{:d}".format(dgtf[k]["exons"][i][0], dgtf[k]["exons"][i][1]-1)
				
				if aid not in acceptors:
					acceptors[aid] = {"gene_id":[], "transcript_id":[], "gene_name":[], "strand":[]}
				
				acceptors[aid]["gene_id"].append(dgtf[k]["info"]["gene_id"])
				acceptors[aid]["transcript_id"].append(dgtf[k]["info"]["transcript_id"])
				acceptors[aid]["gene_name"].append(dgtf[k]["info"]["gene_name"])
				acceptors[aid]["strand"].append(dgtf[k]["info"]["strand"])
				
				# junction is from the last donor to this acceptor
				jid = "{:s}:{:d}:{:d}".format(dgtf[k]["exons"][i][0], dgtf[k]["exons"][i-1][2]+1, dgtf[k]["exons"][i][1]-1)
				
				if jid not in junctions:
					junctions[jid] = {"gene_id":[], "transcript_id":[], "gene_name":[], "strand":[]}

				junctions[jid]["gene_id"].append(dgtf[k]["info"]["gene_id"])
				junctions[jid]["transcript_id"].append(dgtf[k]["info"]["transcript_id"])
				junctions[jid]["gene_name"].append(dgtf[k]["info"]["gene_name"])
				junctions[jid]["strand"].append(dgtf[k]["info"]["strand"])
				
			
	sys.stderr.write("writing output to stdout\n")
	
	for k in sorted(donors):
		lid = list(set(donors[k]["gene_id"]))
		lid.sort()
		gid = ",".join(lid)
		
		lid = list(set(donors[k]["transcript_id"]))
		lid.sort()
		tid = ",".join(lid)
		
		lid = list(set(donors[k]["gene_name"]))
		lid.sort()
		gname = ",".join(lid)

		lid = list(set(donors[k]["strand"]))
		lid.sort()
		strand = ",".join(lid)
				
		print "\t".join([k, "donor", gid, tid, gname, strand])

	for k in sorted(acceptors):
		lid = list(set(acceptors[k]["gene_id"]))
		lid.sort()
		gid = ",".join(lid)
		
		lid = list(set(acceptors[k]["transcript_id"]))
		lid.sort()
		tid = ",".join(lid)
		
		lid = list(set(acceptors[k]["gene_name"]))
		lid.sort()
		gname = ",".join(lid)

		lid = list(set(acceptors[k]["strand"]))
		lid.sort()
		strand = ",".join(lid)
				
		print "\t".join([k, "acceptor", gid, tid, gname, strand])

	for k in sorted(junctions):
		lid = list(set(junctions[k]["gene_id"]))
		lid.sort()
		gid = ",".join(lid)
		
		lid = list(set(junctions[k]["transcript_id"]))
		lid.sort()
		tid = ",".join(lid)
		
		lid = list(set(junctions[k]["gene_name"]))
		lid.sort()
		gname = ",".join(lid)

		lid = list(set(junctions[k]["strand"]))
		lid.sort()
		strand = ",".join(lid)
				
		print "\t".join([k, "junction", gid, tid, gname, strand])

	return(0)



# 
# parse_gtf_attr
#
def parse_gtf_attr(field):

	attrs = {}

	# split into fields
	l1 = field.split(";")

	# split each field into key and value
	l2 = []
	for i in range(len(l1)):
		l2.append(l1[i].split("\""))


	for i in range(len(l2)):
		if len(l2[i]) > 1:
			attrs[l2[i][0].strip()] = l2[i][1].strip()

	return(attrs)		

def load_gtf(fname):
	
	dgtf = {}
	
	# open gtf
	fin = open(fname, "r")
	# loop through and parse annotation into transcript bundles
	for szl in fin:
		ll = szl.strip().split("\t")
		
		if ll[2] != "exon":
			continue
		
		attr = parse_gtf_attr(ll[8])
		attr["strand"] = ll[6]
		tid = attr["transcript_id"]
		if tid not in dgtf:
			dgtf[tid] = {}
			dgtf[tid]["info"] = attr
			dgtf[tid]["exons"] = []
		
		# append exon
		dgtf[tid]["exons"].append([ll[0], int(ll[3]), int(ll[4])])
		
		# make sure the set of exons are sorted
		if len(dgtf[tid]["exons"]) > 1:
			dgtf[tid]["exons"].sort(key=lambda x: x[1])
		
	return(dgtf)
				

# -----------------------------------------------------------------------------
# entry point
# -----------------------------------------------------------------------------


parser = argparse.ArgumentParser(description="Creates a donor/acceptor annotation based on GTF for use in annotating the psi/theta splicing type metrics.")
parser.add_argument('gtf', type=str, help="GTF file")

args = parser.parse_args()

if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("killed it\n")
