#!/usr/bin/env python
#
# fa-unwrap.py
#
# Shawn Driscoll
# 20121128
#
# This script processes the output of MUSCLE and can generate either
# a table of per base A/C/T/G/- counts or a final consensus sequence.
#

import sys
import numpy as np

fin = sys.stdin
seq = ""
all_seqs = []

for szl in fin:
	szl = szl.strip()

	if szl[0] == ">":
		if len(seq) > 0:
			all_seqs.append(seq)

		seq = ""
	else:
		seq = seq + szl

if len(seq) > 0:
	all_seqs.append(seq)


#
# all sequences have been unwrapped - let's check each position and create a consensus matrix
#
seq_len = len(all_seqs[0])
num_seqs = len(all_seqs)
counts = {"A":0, "C":0, "T":0, "G":0, "-":0, "N":0}
consensus_list = []
consensus_seq = ""

#print "A\tC\tT\tG\t-"

for i in range(seq_len):
	pos_counts = counts.copy()
	for j in range(num_seqs):
		pos_counts[all_seqs[j][i]] += 1

	#if pos_counts['-']/(num_seqs*1.0) < 0.95:
	#	print "\t".join(map(str,[pos_counts['A'],pos_counts['C'],pos_counts['T'],pos_counts['G'],pos_counts['-']]))

	if pos_counts['-']/(num_seqs*1.0) < 0.95:
		max_hits = 0
		consensus_base = ""
		for key in pos_counts.keys():
			if pos_counts[key] > max_hits:
				max_hits = pos_counts[key]
				consensus_base = key

		if consensus_base != '-':
			consensus_seq += consensus_base

# done
print ">consensus_seq"
print consensus_seq



