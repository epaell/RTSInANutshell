#!/usr/bin/env python
#
# script that reads in a TLE file, estimates the power output by the beamformer, and plots the result.
#
# Requires:
#   - PyEphem
#   - Matplotlib
#

import os
from optparse import OptionParser
import sys, getopt, string, re
import datetime
from numpy import *

def checkbp(flagtiles, cc, plot_chan, plot_raw):
	sel_offset = 1
	bp_file = "BandpassCalibration_node%03d.dat" %(cc)
	print("Reading %s" %(bp_file))
	
	# open the the file
	try:
		fid = open(bp_file, 'r')
	except IOError:
		print("Can't open file %s for reading." %(bp_file))
		return

	# read the file into an array
	lines = fid.readlines()
	# close the file
	fid.close()

	# check that something was read in
	if len(lines)==0:
		print("Error reading bandpass file: no lines read from the file.")
		return
	elif len(lines)%8 != 1:
		print("Error reading bandpass file: number of lines should be 1 plus a multiple of 8.")
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

	if plot_raw:
		ch = range(0,N_ch)
	else:
		ch = zeros(N_ch)
		for k in range(0,N_ch):
			ch[k] = freq[k]/0.04

	freq_idx = argsort(freq)

	chan_sel = sel_offset + 2*array(freq_idx)

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

	if plot_chan==1:
		freq=ch

	band_start = freq[0]
	quart_bw = ( freq[-1] - band_start ) / 4.0

	# --------------------------------------------------------------------------------- #

	if sel_offset == 1:
		thresh = 0.35

	tilerms = []
	for ant in range(0, N_ant):
		tile = PX_lsq[ant][0]
		valPX = std(fabs(PX_lsq[ant][chan_sel]))
		valPY = std(fabs(PY_lsq[ant][chan_sel]))
		valQX = std(fabs(QX_lsq[ant][chan_sel]))
		valQY = std(fabs(QY_lsq[ant][chan_sel]))

		tilerms.append([tile-1, max([valPX, valPY, valQX, valQY]), valPX, valPY, valQX, valQY])

	for ant in range(0, N_ant):
		tile = int(PX_lsq[ant][0])
#		figure(1);

		val = max(fabs(PX_lsq[ant][chan_sel]))
		if fabs( val - 1.0 ) >= thresh:
			flagtiles[tile-1] += 1

		val = max(fabs(PY_lsq[ant][chan_sel]))
		if val >= thresh:
			flagtiles[tile-1] += 1

		val = max(fabs(QX_lsq[ant][chan_sel]))
		if val >= thresh:
			flagtiles[tile-1] += 1

		val = max(fabs(QY_lsq[ant][chan_sel]))
		if fabs( val - 1.0 ) >= thresh:
			flagtiles[tile-1] += 1

parser = OptionParser()#create parser for program options
parser.add_option("-o", "--obsid",type='string',action="store", dest="obs",help='list of obsids', default='default')
options, args = parser.parse_args()


if options.obs != 'default':
    os.chdir(os.path.join(os.getcwd(),options.obs))
plot_chan = 0
plot_raw  = 0
flagtiles = zeros((128))
for cc in range(24):
	checkbp(flagtiles, cc + 1, plot_chan, plot_raw)
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
flagged=[]#list of tiles to wirte to flagged_tiles.txt
for n, t in tile:
	print("Flag tile %d (%d bad bandpass entries)" %(t, n))
	if n > 10: #make list of tiles with bad bandpasses to write to flagged_tiles.txt
        	flagged.append(t)

flagged.sort()
#write flagged_tiles.txt
flag_file = open('flagged_tiles.txt','w')
for i in range (0,len(flagged)):
    	flg='%d\n'%(flagged[i])
    	flag_file.write(flg)
flag_file.flush()
flag_file.close
print ('flag file written')

