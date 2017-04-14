#!/usr/bin/env python
from mwa_pipe import *

rts = MWAPipe("image")

obslist = get_local_obs_list("1165925976", "1165925976")
#
rts.weighting = "robust"
rts.robustness = -1.0
rts.field_size = 8.0
rts.ncalibrators = 1
rts.template_base = "rts_"
rts.fscrunch = 1
rts.cat_extent = 20.0
rts.cal_cadence = 64
rts.img_cadence = 96
rts.cal_srcs = 30
rts.do_accumulate = False
#
# Test set for transient
rts.min_baseline = 50.0
rts.make_psf = False
rts.cat_extent = 6.0
rts.npeel = 0
rts.niono_cal = 0
rts.img_srcs = 100
rts.update_cal = False
#
rts.rts_bin = "rts_cpu"
rts.img_cat = None
for obsid in obslist:
	# Image the obsid using in-field calibration performed earlier
	rts.image(obsid, cal_id = None, ra_hrs = 1.6504278, dec_deg = -17.9501, dest_image_path = "uvceti")

