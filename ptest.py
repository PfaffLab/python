#!/usr/bin/env python

import sys
from GenomeJunk import TranscriptomeAnnotation as ta

if __name__ == "__main__":
	tannot = ta.TranscriptomeAnnotation()
	tannot.load_refflat(sys.argv[1])
	sys.exit(1)


