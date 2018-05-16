#!/usr/bin/env python
#
# how to use dill to write python objects to a binary data file
#

import dill


fin = open("bwa-quant.py", "r")
dout = open("pythonData.dill", "wb")

for szl in fin:
	ll = szl.strip().split()
	dill.dump(ll, dout)

fin.close()
dout.close()

# read it back in
with (open("pythonData.dill", "rb")) as din:
	while True:
		try:
			ll = dill.load(din)
			print ll
		except EOFError:
			break


