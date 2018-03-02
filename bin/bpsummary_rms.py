#!/usr/bin/env python
#
# script that reads in a TLE file, estimates the power output by the beamformer, and plots the result.
#
# Requires:
#   - PyEphem
#   - Matplotlib
#

import sys, getopt, string, re
import datetime
from numpy import *

def checkbp(flagtiles, cc):
	bp_file = "BandpassCalibration_node%03d.dat" %(cc)
	print("Reading %s" %(bp_file))
	
	# open the the file
	try:
		fid = open(bp_file, 'r')
	except IOError:
		print("Warning: Can't open file %s for reading." %(bp_file))
		return

	# read the file into an array
	lines = fid.readlines()
	# close the file
	fid.close()

	# check that something was read in
	if len(lines)==0:
		print("Error: No lines read from the file.")
		return
	elif len(lines)%8 != 1:
		print("Error: Number of lines should be 1 plus a multiple of 8.")
		return

	for lineIndex in range(1, len(lines), 8):  # get 0, 8, 16, ...
		if lines[lineIndex].find("nan") != -1:
			print("Warning: bandpass file contains 'nan' solutions")
			return

	# initialise the body list
	PX_lsq = []
	PY_lsq = []
	QX_lsq = []
	QY_lsq = []
	PX_fit = []
	PY_fit = []
	QX_fit = []
	QY_fit = []

	tmp1 = lines[0].split(',')
	tmp2 = []
	for val in tmp1:
		tmp2.append( float(val) )
	N_ch = len(tmp2)

	freq = zeros(N_ch)
	for k in range(0,N_ch):
		freq[k] = tmp2[k]

	ch = zeros(N_ch)
	for k in range(0,N_ch):
		ch[k] = freq[k]/0.04

	freq_idx = argsort(freq)

	chan_sel = 1 + 2*array(freq_idx)

	for lineIndex in range(1, len(lines), 8):  # get 0, 8, 16, ...
		tmp = lines[lineIndex+0].split(",")
		PX_lsq.append([])
		PX_lsq[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			PX_lsq[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+1].split(",")
		PX_fit.append([])
		PX_fit[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			PX_fit[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+2].split(",")
		PY_lsq.append([])
		PY_lsq[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			PY_lsq[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+3].split(",")
		PY_fit.append([])
		PY_fit[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			PY_fit[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+4].split(",")
		QX_lsq.append([])
		QX_lsq[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			QX_lsq[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+5].split(",")
		QX_fit.append([])
		QX_fit[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			QX_fit[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+6].split(",")
		QY_lsq.append([])
		QY_lsq[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			QY_lsq[-1][k] = float(tmp[k])
		tmp = lines[lineIndex+7].split(",")
		QY_fit.append([])
		QY_fit[-1]=zeros(N_ch*2+1)
		for k in range(0,len(tmp)):
			QY_fit[-1][k] = float(tmp[k])

	N_ant = len(PX_lsq)

#	print("Examining Jones matrices for %d frequency channels and %d antennas" % (N_ch, N_ant))

	# --------------------------------------------------------------------------------- #

	band_start = freq[0]
	quart_bw = ( freq[-1] - band_start ) / 4.0

	# --------------------------------------------------------------------------------- #

	tilerms = []
	for ant in range(0, N_ant):
		tile = int(PX_lsq[ant][0])
		valPX = std(fabs(PX_lsq[ant][chan_sel]))
		valPY = std(fabs(PY_lsq[ant][chan_sel]))
		valQX = std(fabs(QX_lsq[ant][chan_sel]))
		valQY = std(fabs(QY_lsq[ant][chan_sel]))
		flagtiles[tile-1] += max([valPX, valPY, valQX, valQY])


flagtiles = zeros((128), dtype=float)
for cc in range(24):
	checkbp(flagtiles, cc + 1)

print("Summary:")
tile = []
#nflag = []
for t in range(len(flagtiles)):
	if flagtiles[t] > 0:
		tile.append((flagtiles[t], t))
#		nflag.append(flagtiles[t])
#		print("  Tile %d: %d flags" %(t, flagtiles[t]))
if len(tile) == 0:
	print("Nothing to flag.")
	sys.exit(0)
tile.sort()
tile.reverse()
for n, t in tile:
	print("Flag tile %d (std: %.2f)" %(t, n))
