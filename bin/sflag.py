#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
import astropy.wcs as wcs
import astropy.io.fits as fits

def stats(obsid, maxrmsfactor):
	workdir = os.getcwd()
	os.chdir(obsid)
	os.system("mkdir -p bad")
	os.system("mv bad/* .")
	flist = glob.glob("2*V.fits")
	flist.sort()

	# Get image shape
	hdu = fits.open(flist[0])
	header = hdu[0].header
	nx = int(header['NAXIS1'])
	ny = int(header['NAXIS2'])
	hdu.close()
	
	win = int(0.25 * min(nx, ny))
	nchan = len(flist)
	rms = []
	freqs = []
	for chan in range(nchan):
		freq = flist[chan].split("_")[2]
		ihdulist = fits.open(flist[chan])
		fmin = np.nanmin(ihdulist[0].data[win:-win, win:-win])
		fmax = np.nanmax(ihdulist[0].data[win:-win, win:-win])
		fstd = np.nanstd(ihdulist[0].data[win:-win, win:-win])
		ihdulist.close()
		print("%s %.3f %.3f %.3f" %(flist[chan], fmin, fmax, fstd))
		freqs.append(freq)
		rms.append(fstd)
#		if fstd > maxrms:
#			os.system("mv *%s*.fits bad" %(freq))
	medrms = np.nanmedian(np.array(rms))
	print("Median RMS = %f" %(medrms))
	for chan in range(nchan):
		if rms[chan] > medrms * maxrmsfactor:
			os.system("mv *%s*.fits bad" %(freqs[chan]))
	os.chdir(workdir)

obsids = glob.glob("1?????????")
obsids.sort()
maxrmsfactor = 1.8
if len(sys.argv) > 2:
	maxrmsfactor = float(sys.argv[2])
if len(sys.argv) > 1:
	findex = int(sys.argv[1])
else:
	sys.exit("usage: sflag.py index <maxrmsfactor>")

obsid = obsids[findex]
print("Processing %s" %(obsid))
stats(obsid, maxrmsfactor)

