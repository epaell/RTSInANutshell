#!/usr/bin/env python
#
# script that reads in a TLE file, estimates the power output by the beamformer, and plots the result.
#
# Requires:
#   - PyEphem
#   - Matplotlib
#

import sys, getopt, string, re
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
import datetime

def plotcc(fig, flagtiles, cc, sel_offset, plot_chan, plot_raw):
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

	print("plotting Jones matrices for %d frequency channels and %d antennas" % (N_ch, N_ant))

	# Step through the channels and print spectra

#	fig1 = figure(figsize=(12,8))

	if sel_offset == 1:
#		figure(1)
		fig.canvas.set_window_title('Amp. of the Jones matrices fits (%s)' %(bp_file)) 
	elif sel_offset == 2:
#		figure(1)
		fig.canvas.set_window_title('Phase of the Jones matrices fits (%s)' %(bp_file)) 

		# unwrap 2pi phase jumps
		for ant in range(0, N_ant):
			last_PX_lsq_phase = PX_lsq[ant][2]
			last_PY_lsq_phase = PY_lsq[ant][2]
			last_QX_lsq_phase = QX_lsq[ant][2]
			last_QY_lsq_phase = QY_lsq[ant][2]
			last_PX_fit_phase = PX_fit[ant][2]
			last_PY_fit_phase = PY_fit[ant][2]
			last_QX_fit_phase = QX_fit[ant][2]
			last_QY_fit_phase = QY_fit[ant][2]
			PX_lsq[ant][2] *= 180/pi
			PY_lsq[ant][2] *= 180/pi
			QX_lsq[ant][2] *= 180/pi
			QY_lsq[ant][2] *= 180/pi
			PX_fit[ant][2] *= 180/pi
			PY_fit[ant][2] *= 180/pi
			QX_fit[ant][2] *= 180/pi
			QY_fit[ant][2] *= 180/pi

			for offset in range(4,2*N_ch+1,2):
				last_PX_lsq_phase = PX_lsq[ant][offset]
				last_PY_lsq_phase = PY_lsq[ant][offset]
				last_QX_lsq_phase = QX_lsq[ant][offset]
				last_QY_lsq_phase = QY_lsq[ant][offset]
				last_PX_fit_phase = PX_fit[ant][offset]
				last_PY_fit_phase = PY_fit[ant][offset]
				last_QX_fit_phase = QX_fit[ant][offset]
				last_QY_fit_phase = QY_fit[ant][offset]

				PX_lsq[ant][offset] *= 180/pi
				PY_lsq[ant][offset] *= 180/pi
				QX_lsq[ant][offset] *= 180/pi
				QY_lsq[ant][offset] *= 180/pi
				PX_fit[ant][offset] *= 180/pi
				PY_fit[ant][offset] *= 180/pi
				QX_fit[ant][offset] *= 180/pi
				QY_fit[ant][offset] *= 180/pi

	# --------------------------------------------------------------------------------- #

	if plot_chan==1:
		freq=ch

	band_start = freq[0]
	quart_bw = ( freq[-1] - band_start ) / 4.0

	rc('font', size=8)
	legendfont = matplotlib.font_manager.FontProperties(size=8)

	# --------------------------------------------------------------------------------- #

	if sel_offset == 1:
		thresh = 0.35
	if sel_offset == 2:
		thresh = 30

	tilerms = []
	for ant in range(0, N_ant):
		tile = PX_lsq[ant][0]
		valPX = np.std(fabs(PX_lsq[ant][chan_sel]))
		valPY = np.std(fabs(PY_lsq[ant][chan_sel]))
		valQX = np.std(fabs(QX_lsq[ant][chan_sel]))
		valQY = np.std(fabs(QY_lsq[ant][chan_sel]))

		tilerms.append([tile-1, np.max([valPX, valPY, valQX, valQY]), valPX, valPY, valQX, valQY])
#		print("%3d: %.4f %.4f %.4f %.4f" %(tile-1, valPX, valPY, valQX, valQY))

#	tilerms=sorted(tilerms, key=lambda a_entry: a_entry[1])
#	for t in range(0, N_ant):
#		print("%03d: %.4f %.4f %.4f %.4f %.4f" %(tilerms[t][0], tilerms[t][1], tilerms[t][2], tilerms[t][3], tilerms[t][4], tilerms[t][5]))

	for ant in range(0, N_ant):
		tile = int(PX_lsq[ant][0])
#		figure(1);

		ax = subplot(2,2,1);
		val = max(fabs(PX_lsq[ant][chan_sel]))
		if fabs( val - 1.0 ) < thresh:
			plot(freq[freq_idx],PX_lsq[ant][chan_sel], '-b.')
		else:
			plot(freq[freq_idx],PX_lsq[ant][chan_sel], '-c.') #, label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
#			print('Possible PP flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
			flagtiles[tile-1] += 1

		ax = subplot(2,2,2);
		val = max(fabs(PY_lsq[ant][chan_sel]))
		if val < thresh or sel_offset == 2:
			plot(freq[freq_idx],PY_lsq[ant][chan_sel], '-b.')
		else:
			plot(freq[freq_idx],PY_lsq[ant][chan_sel], '-c.') #, label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
			if sel_offset == 1:
#				print('Possible PQ flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
				flagtiles[tile-1] += 1

		ax = subplot(2,2,3);
		val = max(fabs(QX_lsq[ant][chan_sel]))
		if val < thresh or sel_offset == 2:
			plot(freq[freq_idx],QX_lsq[ant][chan_sel], '-b.')
		else:
			plot(freq[freq_idx],QX_lsq[ant][chan_sel], '-c.') #, label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
			if sel_offset == 1:
#				print('Possible QP flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
				flagtiles[tile-1] += 1

		ax = subplot(2,2,4);
		val = max(fabs(QY_lsq[ant][chan_sel]))
		if fabs( val - 1.0 ) < thresh:
			plot(freq[freq_idx],QY_lsq[ant][chan_sel], '-b.')
		else:
			plot(freq[freq_idx],QY_lsq[ant][chan_sel], '-c.') #, label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
#			print('Possible QQ flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
			flagtiles[tile-1] += 1

	# --------------------------------------------------------------------------------- #

	for ant in range(0, N_ant):

		tile = PX_lsq[ant][0]

#		figure(1);

		ax = subplot(2,2,1);
		val = max(fabs(PX_lsq[ant][chan_sel]))
		if fabs( val - 1.0 ) < thresh:
			plot(freq[freq_idx],PX_fit[ant][chan_sel], '--r')
		else:
			plot(freq[freq_idx],PX_fit[ant][chan_sel], '-k')

		ax = subplot(2,2,2);
		val = max(fabs(PY_lsq[ant][chan_sel]))
		if val < thresh or sel_offset == 2:
			plot(freq[freq_idx],PY_fit[ant][chan_sel], '--r')
		else:
			plot(freq[freq_idx],PY_fit[ant][chan_sel], '-k')

		ax = subplot(2,2,3);
		val = max(fabs(QX_lsq[ant][chan_sel]))
		if val < thresh or sel_offset == 2:
			plot(freq[freq_idx],QX_fit[ant][chan_sel], '--r')
		else:
			plot(freq[freq_idx],QX_fit[ant][chan_sel], '-k')

		ax = subplot(2,2,4);
		val = max(fabs(QY_lsq[ant][chan_sel]))
		if fabs( val - 1.0 ) < thresh:
			plot(freq[freq_idx],QY_fit[ant][chan_sel], '--r')
		else:
			plot(freq[freq_idx],QY_fit[ant][chan_sel], '-k')

	# --------------------------------------------------------------------------------- #

	ax = subplot(2,2,1);
	legend(prop=legendfont)
	ax.xaxis.set_major_locator(MaxNLocator(4))
	ax.yaxis.set_major_locator(MaxNLocator(5))
	title( "Jones matrix element P <- X" )
	if plot_chan==1:
		ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
	if plot_chan==0:
		xlabel( "MHz" )
	if plot_chan==1:
		xlabel( "dumb chan # (doesn't skip flagged channels)" )
	if sel_offset == 1:
		ylabel( "gain relative to band average" )
	if sel_offset == 2:
		ylabel( "degrees" )
	grid(True)

	ax = subplot(2,2,2);
	#legend(loc=[0,0.85],prop=legendfont)
	ax.xaxis.set_major_locator(MaxNLocator(4))
	ax.yaxis.set_major_locator(MaxNLocator(5))
	title( "Jones matrix element P <- Y" )
	if plot_chan==1:
		ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
	if plot_chan==0:
		xlabel( "MHz" )
	if plot_chan==1:
		xlabel( "dumb chan # (doesn't skip flagged channels)" )
	if sel_offset == 1:
		ylabel( "gain relative to band average" )
	if sel_offset == 2:
		ylabel( "degrees" )
	grid(True)

	ax = subplot(2,2,3);
	ax.xaxis.set_major_locator(MaxNLocator(4))
	ax.yaxis.set_major_locator(MaxNLocator(5))
	title( "Jones matrix element Q <- X" )
	if plot_chan==1:
		ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
	if plot_chan==0:
		xlabel( "MHz" )
	if plot_chan==1:
		xlabel( "dumb chan # (doesn't skip flagged channels)" )
	if sel_offset == 1:
		ylabel( "gain relative to band average" )
	if sel_offset == 2:
		ylabel( "degrees" )
	grid(True)

	ax = subplot(2,2,4);
	#legend(loc=[0,0.85],prop=legendfont)
	legend(prop=legendfont)
	ax.xaxis.set_major_locator(MaxNLocator(4))
	ax.yaxis.set_major_locator(MaxNLocator(5))
	title( "Jones matrix element Q <- Y" )
	if plot_chan==1:
		ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
	if plot_chan==0:
		xlabel( "MHz" )
	if plot_chan==1:
		xlabel( "dumb chan # (doesn't skip flagged channels)" )
	if sel_offset == 1:
		ylabel( "gain relative to band average" )
	if sel_offset == 2:
		ylabel( "degrees" )
	grid(True)

#	show()

plot_chan = 0
plot_raw  = 0
flagtiles = np.zeros((128))
with PdfPages('BandpassSummary.pdf') as pdf:
	sel_offset = 1
	for cc in range(24):
		fig = figure(figsize=(8.5,6.5))
		plotcc(fig, flagtiles, cc + 1, sel_offset, plot_chan, plot_raw)
		pdf.savefig()
#		plt.show()
		close()
	print("Summary:")
	for t in range(len(flagtiles)):
		if flagtiles[t] > 0:
			print("  Tile %d: %d flags" %(t, flagtiles[t]))
	sel_offset = 2
	for cc in range(24):
		fig = figure(figsize=(8.5,6.5))
		plotcc(fig, flagtiles, cc + 1, sel_offset, plot_chan, plot_raw)
		pdf.savefig()
#		plt.show()
		close()

	d = pdf.infodict()
	d['Title'] = 'Bandpass Summary'
	d['Author'] = u'Emil Lenc'
	d['Subject'] = 'Summary of Bandpass'
	d['Keywords'] = 'RTS Bandpass'
	d['CreationDate'] = datetime.datetime.today()
	d['ModDate'] = datetime.datetime.today()
