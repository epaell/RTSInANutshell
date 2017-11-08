#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

headerCards = ["GPSTIME", "EXPOSURE", "FILENAME", "MJD", "DATE-OBS", "LST", "HA", "AZIMUTH", "ALTITUDE", "RA", "DEC", "GRIDNUM", "CREATOR", "PROJECT", "MODE", "RECVRS", "DELAYS", "CENTCHAN", "CHANNELS", "CHANSEL", "SUN-ALT", "FINECHAN", "INTTIME", "NAV_FREQ", "NSCANS", "NINPUTS", "NCHANS", "BANDWDTH", "FREQCENT"]

header = {}
for fname in sys.argv[1:]:
	hdulist=fits.open(fname)
	for headerCard in headerCards:
		header[headerCard] = hdulist[0].header[headerCard]
	nInputs = header["NINPUTS"]
	tileData = hdulist[1].data

	print("Project: %s,    Name: %s" %(header["PROJECT"],header["FILENAME"]))
	print("GPS Time: %d" %(header["GPSTIME"]))
	print("Date: %s" %(header["DATE-OBS"]))
	print("LST = %.6f h" %(header["LST"]/15.0))
	print("RA: %.6f h, Dec: %.5f deg" %(header["RA"] / 15.0, header["DEC"]))
	print("Az: %.3f h, Alt: %.2f deg" %(header["AZIMUTH"] / 15.0, header["ALTITUDE"]))
	print("Grid: %d" %(header["GRIDNUM"]))
	print("NSCANS: %d x %.1f (%ds), NINPUTS: %d" %(header["NSCANS"],header["INTTIME"], header["EXPOSURE"], header["NINPUTS"]))
	print("NCHANS: %d x %d kHz = %.2f MHz (Centre = %d, %.2f)" %(header["NCHANS"], header["FINECHAN"], header["BANDWDTH"], header["CENTCHAN"], header["FREQCENT"]))

	print("Flagging:")
	nFlagged = 0
	pre_rts = -1
	for i in range(nInputs):
		nFlagged += tileData[i]["Flag"]
		rts = int(int(tileData[i]["Input"]) / 2.0) + 1
		if rts == pre_rts:
			continue
		if tileData[i]["Flag"] != 0:
			print("RTS %3d; Tile %4s; Rx=%2s; %8s (Flagged)" %(rts-1, tileData[i]["Tile"], tileData[i]["Rx"], tileData[i]["TileName"]))
		else:
			print("RTS %3d; Tile %4s; Rx=%2s; %8s" %(rts-1, tileData[i]["Tile"], tileData[i]["Rx"], tileData[i]["TileName"]))
		pre_rts = rts
	print("%d tiles flagged" %(nFlagged/2))

	badDipoles = 0
	for i in range(nInputs):
		nbad = np.where(tileData[i]["Delays"] == 32)
		badDipoles += len(nbad[0])
	print("Bad dipoles = %d" %(badDipoles))
