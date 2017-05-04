#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
import astropy.wcs as wcs
import astropy.io.fits as fits

# Get the Stokes V image rms as an estimate of image noise
# Assumes that the specified image list has associated Stokes V images with
# the original image: original_x.fits having counterpart names original_V.fits
def v_weight(flist):
	hdu = fits.open(flist[0])
	header = hdu[0].header
	nx = int(header['NAXIS1'])
	ny = int(header['NAXIS2'])
	win = int(0.25 * min(nx, ny))
	nchan = len(flist)
	rms = []
	for chan in range(nchan):
		fvname = flist[chan]
		fvname.replace(fvname[-6:], "V.fits")
		ihdulist = fits.open(fvname)
		fstd = np.nanstd(ihdulist[0].data[win:-win, win:-win])
		rms.append(fstd)
		ihdulist.close()
	rms2 = np.power(np.array(rms), 2.0)
	srms2 = np.sum(rms2)
	return srms2 / rms2

def collapse(flist, dest):
	hdu = fits.open(flist[0])
	header = hdu[0].header
	nx = int(header['NAXIS1'])
	ny = int(header['NAXIS2'])
	nchan = len(flist)
	idata = np.zeros((nchan, ny, nx), dtype=np.float32)
	w = v_weight(flist)
	for chan in range(nchan):
		ihdulist = fits.open(flist[chan])
		idata[chan, :, :] = ihdulist[0].data
		ihdulist.close()
	mdata = np.average(idata, axis=0, weights = w)
	fits.writeto(dest, mdata, header, clobber=True)
	hdu.close()

flist = sys.argv[2:]
flist.sort()
collapse(flist, sys.argv[1])
