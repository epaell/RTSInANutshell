#!/usr/bin/env python

import glob
import sys

prefix = "2"
if len(sys.argv) > 1:
	prefix = sys.argv[1]

gpsids = glob.glob("1?????????")
gpsids.sort()
for gpsid in gpsids:
	flist = glob.glob("%s/%s*.fits" %(gpsid, prefix))
	print("%s (%d)" %(gpsid, len(flist)))

