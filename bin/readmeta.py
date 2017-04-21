#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

headerCards = ["GPSTIME", "EXPOSURE", "FILENAME", "MJD", "DATE-OBS", "LST", "HA", "AZIMUTH", "ALTITUDE", "RA", "DEC", "RAPHASE", "DECPHASE", "SUN-DIST", "MOONDIST", "JUP-DIST", "GRIDNAME", "GRIDNUM", "CREATOR", "PROJECT", "MODE", "RECVRS", "DELAYS", "CALIBRAT", "CENTCHAN", "CHANNELS", "CHANSEL", "SUN-ALT", "FINECHAN", "INTTIME", "NAV_FREQ", "NSCANS", "NINPUTS", "NCHANS", "BANDWDTH", "FREQCENT"]

header = {}
for fname in sys.argv[1:]:
	hdulist=fits.open(fname)
	for headerCard in headerCards:
		header[headerCard] = hdulist[0].header[headerCard]
#	print(header)
	nInputs = header["NINPUTS"]
	# 1 row per tile input
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

	# 1 per flavour of cable
#	digGains = hdulist[2].data

	# 1 row per tile
	#dtype=(numpy.record, [('Input', '>i2'), ('Antenna', '>i2'), ('Tile', '>i2'), ('Pol', 'S1'), ('Rx', '>i2'), ('Slot', '>i2'), ('Flag', '>i2'), ('Length', 'S14'), ('North', '>f4'), ('East', '>f4'), ('Height', '>f4'), ('Gains', '>i2', (24,)), ('Delays', '>i2', (16,)), ('Flavors', 'S10')])
	#for i in range(256):
	print("Flagging:")
	nFlagged = 0
	for i in range(nInputs):
		nFlagged += tileData[i]["Flag"]
		if tileData[i]["Flag"] != 0:
			print("  Input %s; Ant %s; Tile %s Rx: %s: Flagged" %(tileData[i]["Input"], tileData[i]["Antenna"], tileData[i]["Tile"], tileData[i]["Rx"]))
		else:
			print("  Input %s; Ant %s; Tile %s Rx: %s" %(tileData[i]["Input"], tileData[i]["Antenna"], tileData[i]["Tile"], tileData[i]["Rx"]))
#		print("Delays=",tileData[i]["Delays"])
	print("%d inputs flagged" %(nFlagged))

	tiles = np.zeros((int(nInputs/2)), dtype=np.int)
	tiledata = np.zeros((int(nInputs/2)), dtype=np.int)
	badDipoles = 0
	for i in range(nInputs):
	#	print("Input =", tileData[i]["Input"])
	#	print("Antenna =", tileData[i]["Antenna"])
	#	print("Tile =", tileData[i]["Tile"])
	#	print("Pol =", tileData[i]["Pol"])
	#	print("Rx =", tileData[i]["Rx"])
	#	print("Slot =", tileData[i]["Slot"])
	#	print("Flag =", tileData[i]["Flag"])
	#	print("Length =", tileData[i]["Length"])
	#	print("North =", tileData[i]["North"])
	#	print("East =", tileData[i]["East"])
	#	print("Height =", tileData[i]["Height"])
	#	print("Gains =", tileData[i]["Gains"])
#		print("A=%03d T=%03d f=%d %s" %(tileData[i]["Antenna"], tileData[i]["Tile"], tileData[i]["Flag"], tileData[i]["Delays"]))
		nbad = np.where(tileData[i]["Delays"] == 32)
		tiles[tileData[i]["Antenna"]] += len(nbad[0])
		badDipoles += len(nbad[0])
		tiledata[tileData[i]["Antenna"]] += 1
	#	print("Flavors =", tileData[i]["Flavors"])
	nbad = np.where(tiles > 0)
	print("Bad dipoles = %d" %(badDipoles))
	print("Bad tiles = %d / %d" %(len(nbad[0]), len(tiles)))
#	print(tiledata)
sys.exit(0)

inputs = np.zeros((nInputs), dtype=np.int)
ppds = hdulist[3].data
fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in range(ppds.shape[0]):
	plot, = ax1.plot(range(512), ppds[i]["Power"], marker='+', color="black")
plt.show()
plt.close()

#
print()
i=2
print("Input =", mdt[i]["Input"])
print("Antenna =", mdt[i]["Antenna"])
print("Tile =", mdt[i]["Tile"])
print("Pol =", mdt[i]["Pol"])
print("Rx =", mdt[i]["Rx"])
print("Slot =", mdt[i]["Slot"])
print("Flag =", mdt[i]["Flag"])
print("Length =", mdt[i]["Length"])
print("North =", mdt[i]["North"])
print("East =", mdt[i]["East"])
print("Height =", mdt[i]["Height"])
print("Gains =", mdt[i]["Gains"])
print("Delays =", mdt[i]["Delays"])
print("Flavors =", mdt[i]["Flavors"])

# 1 per flavour of cable
dg = hdulist[2].data
print(dg[0])

# 1 per time stamp
# Time, TileID, Power
ppds = hdulist[3].data
for i in range(ppds.shape[0]):
	print(ppds[i]["Time"], ppds[i]["TileID"], len(ppds[i]["Power"]))
#print(ppds[1]["Time"], ppds[1]["TileID"], len(ppds[1]["Power"]))
