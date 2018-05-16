#!/usr/bin/python
#==============================================================================
# template.py
#
# Shawn Driscoll
# date
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# About 
#==============================================================================

import sys
import argparse
import math
import re
from os.path import isfile, expanduser
from collections import defaultdict
from time import localtime
from Basics import messages as ms

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
	
	

	return 0



#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# adding mutually exclusive options
#mygroup = parser.add_mutually_exclusive_group(required=True)
# adding a group

parser.add_argument('fin', type=str, help="Input file")

args = parser.parse_args()

if __name__ == "__main__":

	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nkilled it\n")
	except Exception, e:
		ms.print_exception()
		sys.exit(1)

