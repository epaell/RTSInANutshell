#!/usr/bin/env python
from mwa_pipe import *

rts = MWAPipe("load")

# Get a list of obsids from the archive for UVCeti observed on the 2016-12-18
obslist = get_obs_list(params={'obsname': 'UVCeti_121', 'mintime': 1165925976, 'maxtime': 1165925976})
for obsid in obslist:
	# Fetch the data from the archive
	rts.fetch_data(obsid)
	# Reflag the observation
	rts.reflag(obsid, 1.3)
