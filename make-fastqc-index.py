#!/usr/bin/python
#==============================================================================
# make-fastqc-index.py
#
# Shawn Driscoll
# 20160927
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Builds an index html file for fastqc reports within the current folder
#==============================================================================

import sys, argparse, math, re, os
from os.path import isfile, expanduser

# from igraph import *
# from subprocess import Popen
# from random import gauss, random, sample
# from scipy.stats import norm
# import numpy as np
# import numpy.random as npr

# R support
# import rpy2.robjects as robjects
# r = robjects.r

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
	rres = 0
	d = {} # dict to collect samples and urls
	sid = ""

	# get list of html files at specified path
	rres = os.popen("find {} -name \"fastqc_report.html\"".format(args.path)).readlines()

	for f in rres:
		f = f.strip()
		tmp = f.split("/")
		sid = re.sub("\_fastqc$", "", tmp[-2])
		d[sid] = f

	fout = open(args.o, "w")
	rres = write_header(fout)
	rres = open_list(fout)

	for sid in sorted(d.keys()):
		write_list_item(fout, make_a(sid, d[sid]))

	rres = close_list(fout)
	rres = write_footer(fout)
	fout.close()

	return 0


#==============================================================================
# defs
#==============================================================================

def write_header(f):
	f.write("<!doctype html>\n")
	f.write("<html lang=\"en\">\n")
	f.write("\t<meta charset=\"utf-8\">\n")
	f.write("\t<title>FASTQC Results Index</title>\n")
	f.write("</head>\n\n")
	f.write("<body>\n")
	return 0

def open_list(f):
	f.write("\t<ul>\n")
	return 0

def close_list(f):
	f.write("\t</ul>\n")
	return 0

def make_a(label, url):
	sout = "<a href=\"{}\">{}</a>".format(url, label)
	return sout

def write_list_item(f, txt):
	sout = "<li>{}</li>".format(txt)
	f.write("\t\t" + sout + "\n")
	return 0

def write_footer(f):
	f.write("\n</body>\n")
	f.write("</html>\n")
	return 0


#==============================================================================
# entry point
#==============================================================================


parser = argparse.ArgumentParser(description="Generates an HTML index to link to fastqc report html files in folders")
parser.add_argument('path', type=str, help="Path to FASTQC output folders")
parser.add_argument('-o', type=str, default="index.html", 
	help="Default index file name [index.html]")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")

