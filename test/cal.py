#!/usr/bin/env python
from mwa_pipe import *

rts = MWAPipe("cal")

rts.field_size = 8.0
rts.update_cal = True
rts.ncalibrators = 1
rts.template_base = "rts_"
rts.fscrunch = 1
rts.cat_extent = 20.0
rts.cal_cadence = 64
rts.img_cadence = 96
rts.npeel = 0
rts.niono_cal = 0
rts.cal_srcs = 30
rts.rts_bin = "rts_gpu"
obslist = get_local_obs_list("1165925976", "1165925976")
for obsid in obslist:
	# Clean out any old files
	rts.clean(obsid)
	# Calibrate the obsid
	rts.calibrate(obsid)
