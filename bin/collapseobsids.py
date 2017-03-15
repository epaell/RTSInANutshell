#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
import astropy.wcs as wcs
import astropy.io.fits as fits

def collapse(flist, dest):
	hdu = fits.open(flist[0])
	header = hdu[0].header
	nx = int(header['NAXIS1'])
	ny = int(header['NAXIS2'])
	nchan = len(flist)
	idata = np.zeros((nchan, ny, nx), dtype=np.float32)
	for chan in range(nchan):
		ihdulist = fits.open(flist[chan])
		idata[chan, :, :] = ihdulist[0].data
		ihdulist.close()
	mdata = np.mean(idata, axis=0)
	fits.writeto(dest, mdata, header, clobber=True)
	hdu.close()

obsids = glob.glob("1?????????")
obsids.sort()
findex = int(sys.argv[1])
prefix = "2"
if len(sys.argv) > 2:
	prefix = sys.argv[2]
obsid = obsids[findex]
for stokes in ["I", "Q", "U", "V"]:
	flist = glob.glob("%s/%s*%s.fits" %(obsid, prefix, stokes))
	print("%s_%s_%s.fits" %(prefix, obsid, stokes))
	collapse(flist, "%s_%s_%s.fits" %(prefix, obsid, stokes))
