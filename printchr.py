#!/usr/bin/env python


for i in range(41):
	print chr(33+i) + "\t\t" + "{:0.2%}".format(i/40.0)


