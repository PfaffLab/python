#!/usr/bin/env python
#==============================================================================
# write-bin.py
#
# Shawn Driscoll
# 20120912
#
# Testing out writing data in binary format to a file
#==============================================================================

import sys
from struct import *
from random import randint

chrom = "chr1"
chrom_out = ""
length = 197195432
chrom_field_length = 64
pad = 0
value = 1209



fp = open("test.bin","wb")

#
# write chromosome name
#
pad = chrom_field_length - len(chrom)
fp.write(chrom)
for i in range(pad):
	fp.write(pack("B",0))

#
# write chomoromse length
fp.write(pack("<L",length))

#
# write a bunch of random integers out
for i in range(100):
	fp.write(pack("<H",randint(0,1000)))


fp.close()

