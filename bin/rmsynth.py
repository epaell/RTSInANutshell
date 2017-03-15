#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
import sys
import astropy.wcs as wcs
import astropy.io.fits as fits
import astropy.constants as const

# Perform RM-synthesis on Stokes Q and U data
#
# dataQ, dataU and freqs - contains the Q/U data at each frequency (in Hz) measured.
# startPhi, dPhi - the starting RM (rad/m^2) and the step size (rad/m^2)
def getFDF(dataQ, dataU, freqs, startPhi, stopPhi, dPhi, dType='float32'):
	# Calculate the RM sampling
	phiArr = np.arange(startPhi, stopPhi, dPhi)

	# Calculate the frequency and lambda sampling
	lamSqArr = np.power(const.c.value / np.array(freqs), 2.0)

	# Calculate the dimensions of the output RM cube
	nPhi = len(phiArr)

	# Initialise the complex Faraday Dispersion Function (FDF)
	FDF = np.ndarray((nPhi), dtype='complex')
	
	# Assume uniform weighting
	wtArr = np.ones(len(lamSqArr), dtype=dType)

	K = 1.0 / np.nansum(wtArr)
	
	# Get the weighted mean of the LambdaSq distribution (B&dB Eqn. 32)
	lam0Sq = K * np.nansum(lamSqArr)

	# Mininize the number of inner-loop operations by calculating the
	# argument of the EXP term in B&dB Eqns. (25) and (36) for the FDF
	a = (-2.0 * 1.0j * phiArr)
	b = (lamSqArr - lam0Sq) 
	arg = np.exp( np.outer(a, b) )

	# Create a weighted complex polarised surface-brightness cube
	# i.e., observed polarised surface brightness, B&dB Eqns. (8) and (14)
	Pobs = (np.array(dataQ) + 1.0j * np.array(dataU))

	# Calculate the Faraday Dispersion Function
	# B&dB Eqns. (25) and (36)
	FDF = K * np.nansum(Pobs * arg, 1)
	return FDF, phiArr

def getFreqs(prefix):
	freqs = []
	files = []
	flist = glob.glob("%s*Q.fits" %(prefix))
	for fname in flist:
		pos = fname.find("MHz")
		if pos == -1:
			continue
		freq = fname[pos-7:pos]
		if freq in freqs:
			continue
		freqs.append(freq)
		files.append(fname)
	return freqs, files

def freqInfo(prefix):
	fstrs, fnames = getFreqs(prefix)
	freqs = []
	for fstr in fstrs:
		freqs.append(float(fstr) * 1.0e6)
	freqs.sort()
	fnames.sort()
	df = []
	for f in range(1,len(freqs)):
		df.append(freqs[f] - freqs[f-1])
	return fnames, np.array(freqs), np.median(df)

def findpeak(freqs, fdf, phi, rmsf, rmsfphi, Gauss, dosub):
	components = np.zeros((len(phi)), np.float32)
	peaks = []
	phis = []
	std = 0.0
	rmsflen = int((len(rmsf) - 1) / 2)
	fdflen = len(phi) + rmsflen
	std = np.std(np.abs(fdf))
	peak1 = np.max(np.abs(fdf))
	pos1 = np.argmax(np.abs(fdf))
	val1 = phi[pos1]
	if dosub == False:
		return peak1, val1
	if rmsflen - pos1 < 0:
		return 0.0, 0.0
	fdf -= rmsf[rmsflen - pos1:fdflen - pos1] * fdf[pos1]
	peaks.append(peak1)
	phis.append(val1)
	components[pos1] += peak1
	fdf += np.convolve(components, Gauss, mode='valid')
	return peak1, val1

#-----------------------------------------------------------------------------#
# Main control function                                                       #
#-----------------------------------------------------------------------------#
def RMprocess(prefix, startPhi, stopPhi, dPhi, doclean, writeCube):
	fnames, freqs, df = freqInfo(prefix)
	nchan = len(freqs)
	print("Channel separation = %f" %(df))
	print("Frequency range: %f-%f" %(np.min(freqs), np.max(freqs)))

	# Get the spatial dimensions to use from the first file
	hdu = fits.open(fnames[0])
	header = hdu[0].header
	nx = int(header['NAXIS1'])
	ny = int(header['NAXIS2'])
	
	# Allocate space for the data
	dataQ = np.zeros((nchan, ny, nx), dtype=np.float32)
	dataU = np.zeros((nchan, ny, nx), dtype=np.float32)
	
	for chan in range(nchan):
#		fnames = glob.glob("fr_%7.3fMHz_Q.fits" %(freqs[chan] / 1.0e6))
#		f = freqs[chan]
		qhdulist = fits.open(fnames[chan])
		dataQ[chan, :, :] = qhdulist[0].data
		qhdulist.close()
		uhdulist = fits.open(fnames[chan].replace("Q", "U"))
		dataU[chan, :, :] = uhdulist[0].data
		uhdulist.close()

	# Generate the array of Phi values
	phiArr = np.arange(startPhi, stopPhi, dPhi)
	sys.stdout.flush()

	# Calculate the frequency and lambda sampling from the available frequencies
	#    freqArr_Hz = fits_make_axis(headQ, 2, f0, fc, nc)
	lamArr_m = const.c.value / np.array(freqs)
	lamSqArr_m2 = np.power(lamArr_m, 2.0)

	# Transpose XYZ into ZXY order (spectrum first)
	# Remember Python ordering of arrays is reversed [zyx]
	# Reorder [zyx] -> [yxz], i.e., [0,1,2] -> [1,2,0]
	dataQ = np.transpose(dataQ, (1,2,0))
	dataU = np.transpose(dataU, (1,2,0))
	
	# Run the RM-Synthesis routine on the data
	# Currently assume 3-dimensions in do_rmsynth function
	# do_rmsynth returns a complex FDF cube in spectral order [yxz]
	print(" Running RM-Synthesis routine ...")
	sys.stdout.flush()
	zeroPol = dataQ - dataQ
#	FDFcube = do_rmsynth(dataQ, zeroPol, lamSqArr_m2, phiArr)
#	FDFcube = do_rmsynth(dataU, zeroPol, lamSqArr_m2, phiArr)
	FDFcube = do_rmsynth(dataQ, dataU, lamSqArr_m2, phiArr)

	if doClean == True:
		rstartPhi = startPhi * 2
		rstopPhi = stopPhi * 2 - dPhi
		RMSF, rmsfphi = getFDF(np.ones((len(freqs))), np.zeros((len(freqs))), freqs, rstartPhi, rstopPhi, dPhi)
		# Create the Gaussian filter for reconstruction
		c = 299792458.0 # Speed of light
		lam2 = (c / freqs) ** 2.0
		lam02 = np.mean(lam2)
		minl2 = np.min(lam2)
		maxl2 = np.max(lam2)
		width = (2.0 * np.sqrt(3.0)) / (maxl2 - minl2)
		Gauss = np.exp((-rmsfphi ** 2.0) / (2.0 * ((width / 2.355) ** 2.0)))
		for y in range(ny):
			for x in range(nx):
				fdf = FDFcube[x,y,:]
				findpeak(freqs, fdf, phiArr, RMSF, rmsfphi, Gauss, True)
				FDFcube[x,y,:] = fdf

	# Transpose the data back into image order [yxz] -> [zyx] ([012] -> [201])
	FDFcube = np.transpose(FDFcube, (2,0,1))

	fdf = np.abs(FDFcube)
	datapeak = np.nanmax(fdf, axis=(0))
	dataval = phiArr[np.nanargmax(fdf, axis=(0))]

	# Make a copy of the Q header and alter Z-axis as Faraday depth
	headRM = header.copy()
	headRM['CTYPE3']='FARADAY DEPTH'
	headRM['CDELT3']= phiArr[1] - phiArr[0]
	headRM['CRPIX3']= 1.0
	headRM['CRVAL3']= phiArr[0]

	# If doing a single component clean then add a "c" prefix to the output files to distinguish between cleaned and dirty outputs.
	fprefix = ""
	if doClean == True:
		fprefix = "c"

	# Save a polarized intensity cube
	absFDFcube = np.abs(FDFcube)
	if writeCube == True:
		# write the RM cube (flux density)
		fits.writeto('%sPI.fits' %(fprefix), absFDFcube, headRM, overwrite=True, output_verify='ignore')
		# write the RM cube (rotation angle)
		fits.writeto('%sPIphi.fits' %(fprefix), np.angle(FDFcube), headRM, overwrite=True, output_verify='ignore')
	# Write the peak polarised flux derived from the RM cube
	fits.writeto('%speak.fits' %(fprefix), datapeak, header, overwrite=True, output_verify='ignore')
	# Write the RM at which the peak polarised flux occurs for each pixel in the RM cube
	fits.writeto('%sval.fits' %(fprefix), dataval, header, overwrite=True, output_verify='ignore')


#-----------------------------------------------------------------------------#
# Perform RM-synthesis on Stokes Q and U cubes                                #
# To-Do: Implement efficient masking.                                         #
# To-Do: Allow RM-synthesis on 1D and 2D arrays.                              #
#-----------------------------------------------------------------------------#
def do_rmsynth(dataQ, dataU, lamSqArr, phiArr, dType='float32'):

	# Parse the weight argument
	wtArr = np.ones(lamSqArr.shape, dtype=dType)

	# Sanity check on Q & U data array sizes
	if not dataQ.shape == dataU.shape:
		sys.exit("do_rmsynth: Stokes Q and U data arrays be the same shape.")

	# Check that the 
	if not dataQ.shape[-1] == lamSqArr.shape[-1]:
		sys.exit("do_rmsynth: The Stokes Q and U arrays must be in spectral order.\n # Stokes = %d, # Lamda = %d." % (dataQ.shape[-1], lamSqArr.shape[-1]))

	# Calculate the dimensions of the output RM cube
	nX = dataQ.shape[1]
	nY = dataQ.shape[0]
	nPhi = phiArr.shape[0]

	# Initialise the complex Faraday Dispersion Function (FDF) cube
	# Remember, python index order is reversed [2,1,0] = [y,x,phy]
	FDFcube = np.ndarray((nY, nX, nPhi), dtype='complex')

	# B&dB equations (24) and (38) give the inverse sum of the weights
	K = 1.0 / np.nansum(wtArr)

	# Get the weighted mean of the LambdaSq distribution (B&dB Eqn. 32)
	lam0Sq = K * np.nansum(lamSqArr)

	# Mininize the number of inner-loop operations by calculating the
	# argument of the EXP term in B&dB Eqns. (25) and (36) for the FDF
	a = (-2.0 * 1.0j * phiArr)
	b = (lamSqArr - lam0Sq) 
	arg = np.exp( np.outer(a, b) )

	# Create a weighted complex polarised surface-brightness cube
	# i.e., observed polarised surface brightness, B&dB Eqns. (8) and (14)
	# Weight-array will broadcast to the spectral dimension
	PobsCube = (dataQ + 1.0j * dataU)

	# Do the synthesis at each pixel of the image
	# Note: There may be numpy function which will speed this up
	#       and remove the need for loop. Would also negate the need to
	#       define an empty FDFcube above.
	for k in range(nY):
		for i in range(nX):
			# Calculate the Faraday Dispersion Function
			# B&dB Eqns. (25) and (36)
			FDFcube[k,i,:] = K * np.nansum(PobsCube[k,i,:] * arg, 1)

	return FDFcube

#prefix = "2"
#startPhi = -60.0
#dPhi     = 0.1
#stopPhi = -startPhi+dPhi
doClean = False
writeCube = False
if len(sys.argv) != 7:
	sys.exit("Usage: %s prefix startPhi endPhi dPhi clean cube\nWhere prefix is used to identify image files; startPhi, endPhi and dPhi define the phi range and resolution; clean is 0 or 1 depending on whether the first peak in the FDF should be cleaned; and cube is 0 or 1 depending on whether the RM cube should be written" %(sys.argv[0]))
prefix = sys.argv[1]
startPhi = float(sys.argv[2])
stopPhi = float(sys.argv[3])
dPhi = float(sys.argv[4])
if sys.argv[5] == "1":
	doClean = True
if sys.argv[6] == "1":
	writeCube = True

# RM-range to be searched and the step size
RMprocess(prefix, startPhi, stopPhi, dPhi, doClean, writeCube)
