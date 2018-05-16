#!/usr/bin/env python
#==============================================================================
# make-gtf-db.py
#
# Shawn Driscoll
# 20120723
#
# Generate a SQLite database file from a GTF file
#==============================================================================

import sys,argparse,sqlite3
import HTSeq as hts

#==============================================================================
# parse arguments
#==============================================================================

parser = argparse.ArgumentParser(description="Generate a SQLite database file from a GTF file.")
parser.add_argument('infile',type=str,help="GTF file")

args = parser.parse_args()

# check for file
try:
	fin = open(args.infile,"r")
	fin.close()
except IOError as e:
	print "Error: unable to find/open input file: " + args.infile
	sys.exit()

#==============================================================================
# variables
#==============================================================================

szFilename = args.infile
dbFilename = szFilename + ".db"
nRows = 0
nRowid = 0
ni = 0

#==============================================================================
# main script
#==============================================================================


# create connection
conn = sqlite3.connect(dbFilename)
c = conn.cursor()

# create tables
c.execute("CREATE TABLE IF NOT EXISTS transcripts (id INTEGER PRIMARY KEY AUTOINCREMENT, tname TEXT, gname TEXT, strand TEXT)")
c.execute("CREATE TABLE IF NOT EXISTS exons (id INTEGER PRIMARY KEY AUTOINCREMENT, tid INTEGER, chrom TEXT, start INTEGER, end INTEGER)")
conn.commit()

# read data from file and load into database
gff = hts.GFF_Reader(args.infile)

sys.stderr.write("\n")

for feature in gff:
	if feature.type == "exon":
		ni += 1
		if ni % 100 == 0:
			sys.stderr.write("\r> rows processed: " + str(ni))
		if ni % 10000 == 0:
			conn.commit()

		szTid = feature.attr['transcript_id']
		szGid = feature.attr['gene_id']

		# is this transcript id already inserted?
		nRows = 0
		for row in c.execute("select rowid from transcripts where tname = '" + szTid + "'"):
			nRowid = row[0]
			nRows += 1

		if nRows == 0:
			# insert this transcript into transcripts table
			result = c.execute("INSERT INTO transcripts (tname,gname,strand) VALUES ('" + szTid + "','" + szGid + "','" + feature.iv.strand + "')")
			nRowid = c.lastrowid

		# insert location data
		c.execute("INSERT INTO exons (tid,chrom,start,end) VALUES (" + str(nRowid) + ",'" + feature.iv.chrom + "'," + str(feature.iv.start) + "," + str(feature.iv.end) + ")")

		#conn.commit()

sys.stderr.write("\r> rows processed: " + str(ni) + "\n")

conn.commit()
c.close()
conn.close()



