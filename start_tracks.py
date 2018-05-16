#!/usr/bin/python
#
# 
# start_tracks.py
# Shawn Driscoll
#
# 20160815
#
# This script generates a start point for a ucsc hub trackDb.txt file from the 
# sample bigwig files in whatever the current folder is.
#

import os, sys, re
import argparse


parser = argparse.ArgumentParser(description="Start UCSC hub tracks from bigwig files", 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bigwig', type=str, nargs="+", help='One or more bigwig files')
parser.add_argument('--bin', type=str, default="all", 
	help="Group samples on this field (assuming the files have the usual name formatting). You can specify multiple fields separated by comma with no spaces. Fields are: date, labber, type, geno, age, sort")

args = parser.parse_args()

group_colors = ["27,158,119", 
	"217,95,2", "117,112,179", "231,41,138", 
	"102,166,30", "230,171,2", "166,118,29", 
	"102,102,102"]

group_colors = ["0,205,205","205,41,144","186,186,186","205,38,38","121,205,205","162,205,90","24,116,205","58,95,205","205,85,85","205,96,144","0,205,0","67,205,128","125,38,205","205,51,51","8,8,8","137,104,205","110,110,110","0,0,205","154,205,50","205,155,29","205,198,115","0,154,205","180,82,205"]

if __name__ == "__main__":
	
	# setup variables
	tpl = "track {TRACK_NAME}\nbigDataUrl {FILE_NAME}\nshortLabel {SHORT_LABEL}\nlongLabel {LONG_LABEL}\nparent {PARENT}\ntype bigWig\ncolor {COLOR}"
	
	group_tpl = "track {GROUP_NAME}\ntype bigWig 0 300000\ncompositeTrack on\nshortLabel {SHORT_LABEL}\nlongLabel {LONG_LABEL}\nviewLimits 0:8\ntransformFunc LOG\nvisibility full\nmaxHeightPixels 150:60:11\nshowSubtrackColorOnUi on\nwindowingFunction mean\nconfigurable on\nautoScale off\nalwaysZero on"
		
	dtracks = {}
	cond = "" 
	condPrint = ""
	dinfo = { "date":0, "labber":0, "type":0, "geno":0, "age":0, "sort":0 }
	
	dcond_colors = {}
	cond_idx = 0
	max_colors = len(group_colors)
	
	# get file list
	flist = args.bigwig
	
	if args.bin != "all":
		bin_on_all = False
		bin_list = args.bin.split(",")
		for k in bin_list:
			if k not in dinfo:
				sys.stderr.write("info field specified for binning is invalid: {}\n".format(k))
				sys.exit(1)
	else:
		bin_list = dinfo.keys()
		
	# loop through it
	for f in flist:
		if re.search("\.bw", f):
			
			# ok this is a bigwig
			tmp = f.split(".")
	
			sid, date, labber, type, geno, age, ssort = tmp[0].split("_")

			dinfo = { "date":date, "labber":labber, "type":type, "geno":geno, "age":age, "sort":ssort }
			
			cond_tmp = [dinfo[k] for k in bin_list]
			cond = ".".join(cond_tmp)
						
			if cond not in dcond_colors:
				dcond_colors[cond] = group_colors[cond_idx]
				cond_idx += 1
				if cond_idx >= max_colors:
					cond_idx = 0
				dtracks[cond] = []
		
			condPrint = cond + "_{}".format(sid)
	
			trackline = re.sub("\{TRACK_NAME\}", sid, tpl)
			trackline = re.sub("\{FILE_NAME\}", f, trackline)
			trackline = re.sub("\{SHORT_LABEL\}", condPrint, trackline)
			trackline = re.sub("\{LONG_LABEL\}", "{} {} {} {} {} {} {}".format(sid, date, labber, type, geno, age, ssort), trackline)
			trackline = re.sub("\{PARENT\}", cond, trackline)
			trackline = re.sub("\{COLOR\}", dcond_colors[cond], trackline)
	
			dtracks[cond].append(trackline)
	
	for cond in sorted(dtracks.keys()):
		
		# create the group header
		gline = re.sub("\{GROUP_NAME\}", cond, group_tpl)
		gline = re.sub("\{SHORT_LABEL\}", cond, gline)
		gline = re.sub("\{LONG_LABEL\}", cond, gline)
		print gline + "\n"
		
		for trackline in dtracks[cond]:
			lout = trackline.split("\n")
			for sz in lout:
				print "\t" + sz
			
			print ""
	
