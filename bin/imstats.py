#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
import astropy.wcs as wcs
import astropy.io.fits as fits

def stats(flist):
	flist.sort()
	hdu = fits.open(flist[0])
	header = hdu[0].header
	nx = int(header['NAXIS1'])
	ny = int(header['NAXIS2'])
	win = int(0.25 * min(nx, ny))
	nchan = len(flist)
	idata = np.zeros((nchan, ny, nx), dtype=np.float32)
	rms = []
	for chan in range(nchan):
		ihdulist = fits.open(flist[chan])
		idata[chan, :, :] = ihdulist[0].data
		ihdulist.close()
		fmin = np.nanmin(idata[chan, win:-win, win:-win])
		fmax = np.nanmax(idata[chan, win:-win, win:-win])
		fstd = np.nanstd(idata[chan, win:-win, win:-win])
		rms.append(fstd)
		print("%s %.3f %.3f %.3f" %(flist[chan], fmin, fmax, fstd))
	print("median std = %.3f" %(np.median(np.array(rms))))
stats(sys.argv[1:])
