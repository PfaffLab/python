#!/usr/bin/env python
#==============================================================================
# timelog.py
#
# Shawn Driscoll
# 20170104
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Command line interface to write time stamps to a file along with a comment
#==============================================================================

from time import time,localtime
import sys, argparse
from os.path import isfile, expanduser

#==============================================================================
# globals
#==============================================================================

HOME = expanduser("~")
DEFDIR = "{}/Desktop/tlog.tsv".format(HOME)

#==============================================================================
# main
#==============================================================================

def main(args):
	
	# set variables
	fname = args.o
	tt = None
	ttime = None
	tdate = None
	sout = None
	comment = " ".join(args.message)

	# open output file
	fout = open(fname, "a")

	# get time and format it into a date string and a time string
	tt = localtime()
	tdate = "{}/{}/{}".format(tt.tm_mon, tt.tm_mday, tt.tm_year)
	ttime = "{:02d}:{:02d}:{:02d}".format(tt.tm_hour, tt.tm_min, tt.tm_sec)

	# format output line
	sout = "{}\t{}\t{}\n".format(tdate, ttime, comment)
	# write it to log file
	fout.write(sout)
	# close log file
	fout.close()

	# print logged line out to the command window
	sys.stderr.write("\nLOGGED:\n{}\n".format(sout))

	return 0


#==============================================================================
# entry point
#==============================================================================

# -- parse options and arguments

parser = argparse.ArgumentParser(description="Enter a message and log time and message to a file.")
parser.add_argument('message', type=str, nargs="+", metavar="msg", 
	help="Comment to log with date and time.")
parser.add_argument('-o', type=str, default=DEFDIR, 
	help="Target log file name to append logged time and message.")

# parse
args = parser.parse_args()

# execute main or bail out
if __name__ == "__main__":
	try:
		sys.exit(main(args))
	except KeyboardInterrupt:
		sys.stderr.write("\nExecution terminated.\n")



