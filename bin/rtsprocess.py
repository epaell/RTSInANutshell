#!/usr/bin/env python
import os
import sys
import glob
import logging
import time
import ephem
import numpy as np
import datetime
import pytz
from astropy.io import fits
from mwapy import ephem_utils

CAT_PATH = "/home/elenc/catalog/"
SRC_CAT = "rtscat.txt"
COARSE_CHANNELS = 24

def readMetafits(path):
	header = {}
	radeg = 0.0
	decdeg = 0.0
	logger.info("Reading metafits information for %s" %(path))
	fname = glob.glob("%s/*_metafits_ppds.fits" %(path))[0]
	hdulist=fits.open(fname)
	header["FIELDNAME"] = hdulist[0].header["FILENAME"]
	header["N_SCANS"] = hdulist[0].header["NSCANS"]
	header["N_INPUTS"] = hdulist[0].header["NINPUTS"]
	header["N_CHANS"] = hdulist[0].header["NCHANS"]

	header["INT_TIME"] = hdulist[0].header["INTTIME"]
	header["FREQCENT"] = hdulist[0].header["FREQCENT"]
	header["BANDWIDTH"] = hdulist[0].header["BANDWDTH"]
	header["RA_HRS"] = hdulist[0].header["RA"] / 15.0
	header["LST_HRS"] = hdulist[0].header["LST"] / 15.0
	header["HA_HRS"] = header["LST_HRS"] - header["RA_HRS"]
	header["DEC_DEGS"] = hdulist[0].header["DEC"]
	return header

def getMWAFchans(f):
	hdu = fits.open(f)
	nchan = int(hdu[0].header['NCHANS'])
	hdu.close()
	return nchan

def readMWAFData(obsid):
	# Read flagging information for obsid
	os.system("rm -fr *.mwaf")
	zip_file = "%s_flags.zip" %(obsid)
	if(os.path.isfile(zip_file)):
		os.system("unzip -o -j %s" %(zip_file))
	else:
		logger.warning("Cotter flag zip file for obsid %s does not exist. WARNING: not using Cotter flags." %(obsid))
		return
	mwaf_files = glob.glob("%s_*.mwaf" %(obsid))
	if(len(mwaf_files) != 24):
		logger.warning("Couldn't find 24 MWAF files. WARNING: not using Cotter flags.")
		return
	nchan = getMWAFchans(mwaf_files[0])
	nflagged = np.zeros((COARSE_CHANNELS, nchan), dtype = np.int)
	
	for cc in range(COARSE_CHANNELS):
		infile = "%s_%02d.mwaf" %(obsid, cc + 1)
		hdu = fits.open(infile)

		n_reflag = 0
		flags = hdu[1].data["FLAGS"]
		rows = hdu[1].header["NAXIS2"]
		header = hdu[0].header
		nchan = int(header['NCHANS'])
		nant = int(header['NANTENNA'])
		nbl = int(nant * (nant + 1) / 2)
		ntime = int(rows / nbl)
		dr = flags.reshape(ntime, nbl, nchan)
		nflagged[cc,:] = np.sum(dr, (0, 1))
		hdu.close()
	return nflagged, nchan

def updateMWAFData(obsid, flagged, threshold):
	mwaf_files = glob.glob("%s_*.mwaf" %(obsid))
	if(len(mwaf_files) != 24):
		sys.exit(-1)
	nchan = getMWAFchans(mwaf_files[0])

	base = np.median(flagged)
	bad = np.where(flagged > threshold*base)
	print("Flagging %d channels out of %d" %(len(bad[1]), nchan * COARSE_CHANNELS))
	
	for cc in range(COARSE_CHANNELS):
		infile = "%s_%02d.mwaf" %(obsid, cc + 1)
		logger.info("Updating %s:" %(infile))
		hdu = fits.open(infile, mode="update")

		n_reflag = 0
		for i in range(nchan):
			if flagged[cc][i] > base * threshold:
				logger.info(i)
				hdu[1].header['REFLG_%02d' % n_reflag] = i
				n_reflag += 1

		hdu.close()

def removePols(pols = ["XX", "YY", "XYim", "XYre"]):
	# Try to get a list of all of the instrumental XX images
	flist = glob.glob("2*_*%s.fits" %(pols[0]))
	while len(flist) > 0:
		# Get the frequency of the first item in the list
		freq0 = flist[0].split("_")[2]
		# Find all time stamps available for that frequency (need to minimize the number of files being removed at any time)
		flist = glob.glob("2*%s*%s.fits" %(freq0, pols[0]))
		for fname in flist:
			tindex = fname.split("_")[-1]
			os.system("rm 2*%s" %(tindex))
			os.system("rm 2*%s" %(tindex.replace(pols[0], pols[1])))
			os.system("rm 2*%s" %(tindex.replace(pols[0], pols[2])))
			os.system("rm 2*%s" %(tindex.replace(pols[0], pols[3])))
	# Remove any files
	os.system("rm 2*%s.fits" %(pols[1]))
	os.system("rm 2*%s.fits" %(pols[2]))
	os.system("rm 2*%s.fits" %(pols[3]))
	return

def moveImages(pols, dest): #"%s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
	# Try to get a list of all of the instrumental XX images
	flist = glob.glob("2*_*%s.fits" %(pols[0]))
	while len(flist) > 0:
		# Get the frequency of the first item in the list
		freq0 = flist[0].split("_")[2]
		# Find all time stamps available for that frequency (need to minimize the number of files being removed at any time)
		flist = glob.glob("2*%s*%s.fits" %(freq0, pols[0]))
		for fname in flist:
			tindex = fname.split("_")[-1]
			os.system("mv 2*%s %s" %(tindex, dest))
			os.system("mv 2*%s %s" %(tindex.replace(pols[0], pols[1]), dest))
			os.system("mv 2*%s %s" %(tindex.replace(pols[0], pols[2]), dest))
			os.system("mv 2*%s %s" %(tindex.replace(pols[0], pols[3]), dest))
	# Remove any files
	os.system("mv 2*%s.fits %s" %(pols[1], dest))
	os.system("mv 2*%s.fits %s" %(pols[2], dest))
	os.system("mv 2*%s.fits %s" %(pols[3], dest))

def GenerateRTSInputFiles(obsid, params, dataPath, basename, template, calSrc, rahrs, decdeg, mode = None):
	array = params["array"]
	if not (array == '32T' or array == '128T'):
		logger.error("Array Parameter should be either '32T' or '128T'")
		return
	
	# Set the time.
	logger.info("obsid=%s" %(obsid))
	obstime = ephem_utils.MWATime(gpstime = float(obsid))
	lsthrs = float(obstime.LST) / 15.0
	jd = obstime.MJD + 2400000.5
	logger.info("lst=%.4f hr" %(lsthrs))
	logger.info("jd=%.5f" %(jd))
	header = readMetafits(dataPath)
	bandwidth = header["BANDWIDTH"] / 24
	if header["HA_HRS"] > 12.0:
		header["HA_HRS"] -= 24.0

	# How many GPU boxes are there?
	GPUFiles = glob.glob(dataPath + '/*gpubox*00.fits')  

	if(len(GPUFiles) > 0):
		bandList = [int(filename[-10:-8]) for filename in GPUFiles]
	else:
		uvfileList = glob.glob(dataPath + '/*uvfits')
		bandList = [int(filename[-9:-7]) for filename in uvfileList]   
	bandList.sort()
	if params["maxsubbands"] != "-1":
		bandList = bandList[:int(params["maxsubbands"])]	# Limit the number of subbands to MaxSubBands
	subbandString = ''
	for band in bandList:
		subbandString = subbandString + str(band) + ','
	subbandString = subbandString[:-1]

	corrDumpTime = header["INT_TIME"]

	# Generate a channel flag file to remove centre and edge channels
	templateFile = open(template)
	setSubbands = 0

	# Write a new .in file
	outName = "rts_%s.in" %(basename)
	outFile = open(outName, "w+")
	
	corrDumpsPerCadence = -1.0
	integrationBins = 0.0
	if mode == "calibrate":
		outFile.write("generateDIjones=1\n")
		outFile.write("useStoredCalibrationFiles=0\n")
		imagingCadence = int(params["calcadence"])
		integrationBins = 6 - int(np.log2(corrDumpTime))
	else:
		imagingCadence = int(params["imgcadence"])
		
		integrationBins = 4 - int(np.log2(corrDumpTime))
	corrDumpsPerCadence = (int)(imagingCadence / corrDumpTime)
	if corrDumpsPerCadence == 1:
		integrationBins = 1
	elif corrDumpsPerCadence == 2:
		integrationBins = 2
	elif corrDumpsPerCadence < 8:
		integrationBins = 3
	if mode == "uv":
		outFile.write("generateDIjones=0\n")
		outFile.write("useStoredCalibrationFiles=1\n")
	if mode == "image":
		outFile.write("generateDIjones=0\n")
		outFile.write("useStoredCalibrationFiles=1\n")
		if params["weighting"] == "uniform":
			outFile.write("imgUniformWeighting=1\n")
		elif params["weighting"] == "robust":
			outFile.write("imgRobustWeighting=1\n")
			outFile.write("imgRobustnessParameter=%s\n" %(params["robustness"]))
		if params["doaccumulate"] == "true":
			outFile.write("UseAccumUVsampling=1\n")
		if params["regrid"] == "true":
			outFile.write("MakeImage=\n")
			outFile.write("FieldOfViewDegrees=%s\n" %(params["fieldsize"]))
			outFile.write("RegridMethod=%s\n" %(params["regridmethod"]))
			outFile.write("DoRegriddingWithProjection=%s\n" %(params["nside"]))
		else:
			outFile.write("FieldOfViewDegrees=%s\n" %(params["fieldsize"]))
			if params["makestokes"] == "true":
				outFile.write("MakeStokesSnapshots=\n")
			else:
				outFile.write("MakeImage=\n")
		if params["makepsf"] == "true":
			outFile.write("ImagePSF=1\n")
	if mode == "accumulate":
                outFile.write("generateDIjones=0\n")
                outFile.write("useStoredCalibrationFiles=1\n")
		outFile.write("DoUVsamplingAccum=1\n")
		outFile.write("FieldOfViewDegrees=%s\n" %(params["fieldsize"]))
		outFile.write("MakeStokesSnapshots=\n")
	if mode == "image" or mode == "accumulate":
		if params["shorttaper"] != "none":
			outFile.write("imgShortBaselineTaper=%s\n" %(params["shorttaper"]))
		if params["longtaper"] != "none":
			outFile.write("imgLongBaselineTaper=%s\n" %(params["longtaper"]))
		if params["minbaseline"] != "none":
			outFile.write("imgBaselineMin=%s\n" %(params["minbaseline"]))
		if params["maxbaseline"] != "none":
			outFile.write("imgBaselineMax=%s\n" %(params["maxbaseline"]))

	logger.info("Imaging cadence = %f" %(imagingCadence))
	scanTime = header["N_SCANS"] * header["INT_TIME"] - header["INT_TIME"] * 4
	logger.info("scanTime = %d" %(scanTime))
	nIterations = int(scanTime / imagingCadence)
	logger.info("nIterations = %d" %(nIterations))

	templateFile = open(template)
	for line in templateFile:
		newLine = line
		if (params["startat"] != "0") and (line.find('StartProcessingAt') == 0):
			newLine = "StartProcessingAt=%s\n" %(params["startat"])
		if (params["startat"] != "0") and (line.find('StartIntegrationAt') == 0):
			newLine = "StartIntegrationAt=%s\n" %(params["startat"])
		if line.startswith('BaseFilename'):
			newLine = "BaseFilename=%s/*_gpubox\n" %(dataPath)
			if params["useflag"] == "true":
				newLine += "ImportCotterFlags=1\nImportCotterBasename=%s/%s\n\n" %(dataPath, obsid)
			if params["usemeta"] == "true":
				newLine += "ReadMetafitsFile=1\nMetafitsFilename=%s/%s\n" %(dataPath, obsid)
		if line.startswith('CorrDumpTime'):
			newLine = line.replace(line[len('CorrDumpTime='):],"%.1f" %(header["INT_TIME"]) + '\n')
		if line.startswith('FscrunchChan'):
			newLine = "FscrunchChan=%s\n" %(params["fscrunch"])
		if line.startswith('ChannelBandwidth'):
			newLine = line.replace(line[len('ChannelBandwidth='):],"%.3f" %(header["BANDWIDTH"]/header["N_CHANS"]) + '\n')
		if line.startswith('NumberOfChannels'):
			newLine = line.replace(line[len('NumberOfChannels='):],str(header["N_CHANS"]/24) + '\n')
		if line.startswith('ObservationFrequencyBase'):
			newLine = line.replace(line[len('ObservationFrequencyBase='):],str(header["FREQCENT"] - header["BANDWIDTH"] / 2.0 - 0.02) + '\n')
		if line.startswith('ObservationPointCentreHA'):
			newLine = line.replace(line[len('ObservationPointCentreHA='):],str(header["HA_HRS"]) + '\n')
		if line.startswith('ObservationPointCentreDec'):
			newLine = line.replace(line[len('ObservationPointCentreDec='):],str(header["DEC_DEGS"]) + '\n')
		if ('SourceCatalogueFile' in line) and (calSrc == "self"):
			logger.info("Setting catalogue file for self calibration")
			newLine = "SourceCatalogueFile=catalogue.txt\n"
		if params["adjustcalparams"] == "true":
			if 'NumberOfSourcesToPeel' in line:
					newLine = "NumberOfSourcesToPeel=%s\n" %(params["npeel"])
			if 'NumberOfIonoCalibrators' in line:
					newLine = "NumberOfIonoCalibrators=%s\n" %(params["nionocal"])
			if 'UpdateCalibratorAmplitudes' in line:
					newLine = "UpdateCalibratorAmplitudes=%s\n" %(params["updatecal"])
			if 'PrimaryCalibrator' in line:
				if calSrc == None or calSrc == "self":
					continue
				else:
					newLine = "PrimaryCalibrator=%s" %(calSrc)

		if line.startswith('ObservationImageCentreRA'):
			if rahrs == None:
				newLine = line.replace(line[len('ObservationImageCentreRA='):],str(header["RA_HRS"]) + '\n')
			else:
				if rahrs == "zenith":
					newLine = "ObservationImageCentreRA=%s\n" %(lsthrs)
				else:
					newLine = "ObservationImageCentreRA=%s\n" %(rahrs)
		if line.startswith('ObservationImageCentreDec'):
			if decdeg == None:
				newLine = line.replace(line[len('ObservationImageCentreDec='):],str(header["DEC_DEGS"]) + '\n')
			else:
				if decdeg == "zenith":
					newLine = "ObservationImageCentreDec=%s\n" %(-26.0 - 42.0/60.0 - 11.95/3600.0)
				else:
					newLine = "ObservationImageCentreDec=%s\n" %(decdeg)
		if line.startswith('PrimaryCalibrator'):
			if calSrc == None:
				continue
			else:
				newLine = "PrimaryCalibrator=%s" %(calSrc)
		if line.startswith('CorrDumpsPerCadence'):
			newLine = "CorrDumpsPerCadence=%d\n" %(corrDumpsPerCadence)
		if line.startswith('NumberOfIntegrationBins'):
			newLine = "NumberOfIntegrationBins=%d\n" %(integrationBins)
		if line.startswith('NumberOfIterations'):
			if params["endat"] != "0":
				newLine = "NumberOfIterations=%s\n" %(params["endat"])
			else:
				if(array=='32T'):
					# CorrDumpsPerCadence is currently hard set to 4
					newLine = line.replace(line[len('NumberOfIterations='):],str(header["N_SCANS"] / 4) + '\n')
				else:
					newLine = line.replace(line[len('NumberOfIterations='):],str(nIterations) + '\n')
		if line.startswith('SubBandIDs='):
			newLine = line.replace(line[len('SubBandIDs='):],subbandString)
			setSubbands = 1
		if line.startswith('ObservationTimeBase='):
			newLine = "ObservationTimeBase=%.5f\n" %(jd)
#			newLine = line.replace(line[len('ObservationTimeBase='):], "%.5f" %(jd))	
		outFile.write(newLine)
	if(array=='128T'):
		if(setSubbands == 0):
			outFile.write('SubBandIDs='+subbandString+'\n')
	outFile.close()
	templateFile.close()
	return len(bandList)

def extractArgs(str):
	args = []
	if "(" in str:
		data = str.split("(")
		args.append(data[0].strip())
		if "," in data[1]:
			params = data[1][:-1].split(",")
			for p in params:
				args.append(p.strip())
		else:
			args.append(data[-1][:-1].strip())
	else:
		args.append(str.strip())
	return args

class Clean:
	def __init__(self, targetId):
		self.targetId = targetId

	def execute(self, params):
		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))

		# Remove old files
		logger.info("Cleaning up %s/%s" %(params["workdir"], self.targetId))
		os.system("rm -fr *.log 2*.fits *.btr uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits *.dat cat*.txt rts.in")

class Calibrate:
	def __init__(self, targetId):
		self.targetId = targetId

	def execute(self, params):
		logger.info("Calibrate %s" %(self.targetId))

		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))

		# Remove old files
		logger.info("Removing temporary files in %s/%s" %(params["workdir"], self.targetId))
		os.system("rm -fr *.log 2*.fits *.btr uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits *.dat")

		# Use selfcal of field to do calibration
		logger.info("Generating initial calibration with %s" %(self.targetId))
		if params["calcat"] == "regenerate":
			if os.path.exists("%s_metafits_ppds.fits" %(self.targetId)) == True:
				logger.info("Using %s_metafits_ppds.fits as input" %(self.targetId))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -m %s_metafits_ppds.fits -c %s -s %s/%s -n %s" %(self.targetId, params["catextent"], params["catpath"], params["srccat"], params["calsrcs"]))
			else:
				logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(self.targetId))
			os.system("mv *patch*.txt catalogue.txt")
		else:
			# Re-use a catalogue generated for the specified obsid
			os.system("cp %s/%s/catalogue.txt ." %(params["workdir"], params["calcat"]))

		# Generate an RTS input file based on the template
		logger.info("Generating input files %s" %(self.targetId))
		ncoarse = GenerateRTSInputFiles(self.targetId, params, "%s/%s" %(params["workdir"], self.targetId), "cal", "%s/%scal.in" %(params["workdir"], params["templatebase"]), "self", None, None, "calibrate")
		# Copy over the default flag files
#				os.system("cp %s/fl*.txt ." %(params["workdir"]))

		# Do the calibration
		logger.info("Performing calibration for %s" %(self.targetId))
		os.system("aprun -N 1 -n %d %s rts_cal.in" %(ncoarse + 1, params["rtsbin"]))
			
		logger.info("Finished calibrating %s" %(self.targetId))

class Image:
	def __init__(self, targetId, calId, destImagePath, rahrs, decdeg):
		self.targetId = targetId
		self.calId = calId
		self.destImagePath = destImagePath
		self.rahrs = rahrs
		self.decdeg = decdeg

	def execute(self, params):
		logger.info("Imaging %s" %(self.targetId))
		logger.info("Using calibration from %s" %(self.calId))
		if params["weighting"] == "robust":
			logger.info("Weighting = %s (robustness=%s)" %(params["weighting"], params["robustness"]))
		else:
			logger.info("Weighting = %s" %(params["weighting"]))

		if self.rahrs != None and self.decdeg != None:
			logger.info("Using (%s, %s) as imaging centre" %(self.rahrs, self.decdeg))
		# Move to the target directory
		logger.info("Cleaning up %s/%s" %(params["workdir"], self.targetId))
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		# Remove old files
		removePols()
		removePols(["I", "Q", "U", "V"])
		os.system("rm -fr *.log *.btr pee*.txt rest*.txt cor*.fits int*.fits")
		# Move to final cal to see if there is anything there
		if self.calId == "self":
			os.chdir("%s/%s" %(params["workdir"], self.targetId))
			# Use calibration already in target
		else:
			# Use calibration from specified obsid
			os.chdir("%s/%s" %(params["workdir"], self.calId))
			if self.calId != self.targetId:
				# Copy the calibration solutions to the target
				logger.info("Copying calibration files from %s to %s" %(self.calId, self.targetId))
				# copy across the calibration data ... unless we are imaging the calibrator.
				os.system("cp *.dat %s/%s" %(params["workdir"], self.targetId))

		# Check if a calibration exists
		if len(glob.glob("*.dat")) == 0:
			logger.warning("No calibration data found!")
			return

		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		# Generate the catalogue for imaging.
		if params["imgcat"] == "regenerate":
			if os.path.exists("%s_metafits_ppds.fits" %(self.targetId)) == True:
				logger.info("Using %s_metafits_ppds.fits as input" %(self.targetId))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s_metafits_ppds.fits -c %s -s %s/%s -n %s" %(self.targetId, params["catextent"], params["catpath"], params["srccat"], params["imgsrcs"]))
			else:
				logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(self.targetId))
			os.system("mv *peel*.txt catalogue_imaging.txt")
		else:
			logger.info("Using catalogue from %s" %(params["imgcat"]))
			os.system("cp %s/%s/catalogue_imaging.txt ." %(params["workdir"], params["imgcat"]))
			
		# Generate an RTS input file based on the template
		logger.info("Generating RTS imaging input files for %s" %(self.targetId))
		ncoarse = GenerateRTSInputFiles(self.targetId, params, "%s/%s" %(params["workdir"], self.targetId), "img", "%s/%simg.in" %(params["workdir"], params["templatebase"]), None, self.rahrs, self.decdeg, "image")
		# Do the imaging
		logger.info("Imaging %s" %(self.targetId))
		os.system("aprun -N 1 -n %d %s rts_img.in" %(ncoarse + 1, params["rtsbin"]))
		if params["removeinstpols"] == "true":
			removePols()
		logger.info("Moving image files to destination path: %s" %(self.destImagePath))
		os.system("mkdir -p %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		moveImages(["XX", "YY", "XYim", "XYre"], "%s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		moveImages(["I", "Q", "U", "V"], "%s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
#		os.system("mv *MHz*.fits %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("mv peeled_sources_*.txt %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("mv restore_*.txt %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("mv rts_*.log %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("cp *.dat %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("cp *.in %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		logger.info("Finished imaging %s" %(self.targetId))

class Accumulate:
	def __init__(self, targetId, calId, destImagePath, rahrs, decdeg):
		self.targetId = targetId
		self.calId = calId
		self.destImagePath = destImagePath
		self.rahrs = rahrs
		self.decdeg = decdeg

	def execute(self, params):
		logger.info("Accumulating weights for %s" %(self.targetId))
		logger.info("Using calibration from %s" %(self.calId))

		if self.rahrs != None and self.decdeg != None:
			logger.info("Using (%s, %s) as imaging centre" %(self.rahrs, self.decdeg))
		# Move to the target directory
		logger.info("Cleaning up %s/%s" %(params["workdir"], self.targetId))
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		# Remove old files
		os.system("rm -fr *.log 2*.fits *.btr uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits")
		# Move to final cal to see if there is anything there
		if self.calId == "self":
			os.chdir("%s/%s" %(params["workdir"], self.targetId))
			# Use calibration already in target
		else:
			# Use calibration from specified obsid
			os.chdir("%s/%s" %(params["workdir"], self.calId))
			if self.calId != self.targetId:
				# Copy the calibration solutions to the target
				logger.info("Copying calibration files from %s to %s" %(self.calId, self.targetId))
				# copy across the calibration data ... unless we are imaging the calibrator.
				os.system("cp *.dat %s/%s" %(params["workdir"], self.targetId))

		# Check if a calibration exists
		if len(glob.glob("*.dat")) == 0:
			logger.warning("No calibration data found!")
			return

		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		# Generate the catalogue for imaging.
		if params["imgcat"] == "regenerate":
			if os.path.exists("%s_metafits_ppds.fits" %(self.targetId)) == True:
				logger.info("Using %s_metafits_ppds.fits as input" %(self.targetId))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s_metafits_ppds.fits -c %s -s %s/%s -n %s" %(self.targetId, params["catextent"], params["catpath"], params["srccat"], params["imgsrcs"]))
			else:
				logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(self.targetId))
			os.system("mv *peel*.txt catalogue_imaging.txt")
		else:
			logger.info("Using catalogue from %s" %(params["imgcat"]))
			os.system("cp %s/%s/catalogue_imaging.txt ." %(params["workdir"], params["imgcat"]))
		
		# Generate an RTS input file based on the template
		logger.info("Generating RTS input files to accumulate weights for %s" %(self.targetId))
		ncoarse = GenerateRTSInputFiles(self.targetId, params, "%s/%s" %(params["workdir"], self.targetId), "acc", "%s/%simg.in" %(params["workdir"], params["templatebase"]), None, self.rahrs, self.decdeg, "accumulate")
		# Generate the weights file
		logger.info("Accumulate weights for %s" %(self.targetId))
		os.system("aprun -N 1 -n %d %s rts_acc.in" %(ncoarse + 1, params["rtsbin"]))

		os.system("mv rts_*.log %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("cp *.dat %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		logger.info("Finished accumulating weights for %s" %(self.targetId))

class UVdump:
	def __init__(self, targetId, calId, destImagePath, rahrs, decdeg):
		self.targetId = targetId
		self.calId = calId
		self.destImagePath = destImagePath
		self.rahrs = rahrs
		self.decdeg = decdeg

	def execute(self, params):
		logger.info("UV Dump for %s" %(self.targetId))
		logger.info("Using calibration from %s" %(self.calId))

		if self.rahrs != None and self.decdeg != None:
			logger.info("Using (%s, %s) as imaging centre" %(self.rahrs, self.decdeg))
		# Move to the target directory
		logger.info("Cleaning up %s/%s" %(params["workdir"], self.targetId))
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		# Remove old files
		os.system("rm -fr *.log 2*.fits *.btr uvdump*.fits uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits")
		# Move to final cal to see if there is anything there
		if self.calId == "self":
			os.chdir("%s/%s" %(params["workdir"], self.targetId))
			# Use calibration already in target
		else:
			# Use calibration from specified obsid
			os.chdir("%s/%s" %(params["workdir"], self.calId))
			if self.calId != self.targetId:
				# Copy the calibration solutions to the target
				logger.info("Copying calibration files from %s to %s" %(self.calId, self.targetId))
				# copy across the calibration data ... unless we are imaging the calibrator.
				os.system("cp *.dat %s/%s" %(params["workdir"], self.targetId))

		# Check if a calibration exists
		if len(glob.glob("*.dat")) == 0:
			logger.warning("No calibration data found!")
			return

		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		# Generate the catalogue for imaging.
		if params["imgcat"] == "regenerate":
			if os.path.exists("%s_metafits_ppds.fits" %(self.targetId)) == True:
				logger.info("Using %s_metafits_ppds.fits as input" %(self.targetId))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s_metafits_ppds.fits -c %s -s %s/%s -n %s" %(self.targetId, params["catextent"], params["catpath"], params["srccat"], params["imgsrcs"]))
			else:
				logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(self.targetId))
			os.system("mv *peel*.txt catalogue_imaging.txt")
		else:
			logger.info("Using catalogue from %s" %(params["imgcat"]))
			os.system("cp %s/%s/catalogue_imaging.txt ." %(params["workdir"], params["imgcat"]))
		
		# Generate an RTS input file based on the template
		logger.info("Generating RTS input files to accumulate weights for %s" %(self.targetId))
		ncoarse = GenerateRTSInputFiles(self.targetId, params, "%s/%s" %(params["workdir"], self.targetId), "uv", "%s/%suv.in" %(params["workdir"], params["templatebase"]), None, self.rahrs, self.decdeg, "uv")
		# Generate the weights file
		logger.info("UV Dump for %s" %(self.targetId))
		os.system("aprun -N 1 -n %d %s rts_uv.in" %(ncoarse + 1, params["rtsbin"]))

                os.system("mkdir -p %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("mv rts_*.log %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("cp *.dat %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		os.system("mv uvdump*.uvfits %s/%s/%s" %(params["workdir"], self.destImagePath, self.targetId))
		logger.info("Finished dumping uv for %s" %(self.targetId))

class SetParam:
	def __init__(self, param, value):
		self.param = param
		self.value = value

	def execute(self, params):
		logger.info("Setting parameter %s = %s" %(self.param, self.value))
		params[self.param] = self.value

class Reflag:
	def __init__(self, targetId, threshold):
		self.targetId = targetId
		self.threshold = float(threshold)

	def execute(self, params):
		logger.info("Reflagging %s : threshold=%f" %(self.targetId, self.threshold))
		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		flagged, nchan = readMWAFData(self.targetId)
		updateMWAFData(self.targetId, flagged, self.threshold)
		logger.info("Finished reflagging %s" %(self.targetId))

class PairReflag:
	def __init__(self, targetId, target2, threshold):
		self.targetId = targetId
		self.target2 = target2
		self.threshold = float(threshold)

	def execute(self, params):
		logger.info("Reflagging Pair %s+%s : threshold=%f" %(self.targetId, self.target2, self.threshold))
		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		flagged1, nchan1 = readMWAFData(self.targetId)
		os.chdir("%s/%s" %(params["workdir"], self.target2))
		flagged2, nchan2 = readMWAFData(self.target2)
		if nchan1 != nchan2:
			logger.error("Channels in specified obsids do not match %s=%d; %s=%d" %(obsid1, nchan1, obsid2, nchan2))
			return
		flagged = flagged1 + flagged2
		nchan = nchan1
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		updateMWAFData(self.targetId, flagged, self.threshold)
		os.chdir("%s/%s" %(params["workdir"], self.target2))
		updateMWAFData(self.target2, flagged, self.threshold)
		logger.info("Finished reflagging %s+%s" %(self.targetId, self.target2))

class FetchData:
	def __init__(self, targetId):
		self.targetId = targetId

	def execute(self, params):
		logger.info("Fetch data for %s" %(self.targetId))
		# Move to the target directory
		os.chdir("%s" %(params["workdir"]))
#		os.system("change_db.py curtin")
		os.system("obsdownload.py -o %s" %(self.targetId))
#                os.system("unzip %s_flags.zip" %(self.targetId))
		logger.info("Finished fetching data for %s" %(self.targetId))

class FetchMetadata:
	def __init__(self, targetId):
		self.targetId = targetId

	def execute(self, params):
		logger.info("Fetch metadata for %s" %(self.targetId))
		# Move to the target directory
		os.chdir("%s/%s" %(params["workdir"], self.targetId))
		os.system("wget -O %s_metafits_ppds.fits http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s" %(self.targetId, self.targetId))
	
def readProcessFile(fname):
	items = []
	parameters = {}
	# Set up some default parameters
	catpath = os.environ.get('RTS_CAT_PATH')
	if catpath == None:
		logger.warning("RTS_CAT_PATH not set, using %s as default" %(CAT_PATH))
		catpath = CAT_PATH
	items.append(SetParam("catpath", catpath))
	items.append(SetParam("weighting", "natural"))
	items.append(SetParam("robustness", "+2.0"))
	items.append(SetParam("doaccumulate", "false"))
	items.append(SetParam("workdir", os.getcwd()))
	items.append(SetParam("array", "128T"))
	items.append(SetParam("fieldsize", "2"))
	items.append(SetParam("nside", "2048"))
	items.append(SetParam("maxsubbands", "-1"))
	items.append(SetParam("fscrunch", "4"))

	items.append(SetParam("adjustcalparams", "true"))
	items.append(SetParam("ncalibrators", "1"))
	items.append(SetParam("npeel", "0"))
	items.append(SetParam("nionocal", "0"))
	items.append(SetParam("updatecal", "0"))
	items.append(SetParam("shorttaper", "none"))
	items.append(SetParam("longtaper", "none"))
	items.append(SetParam("minbaseline", "none"))
	items.append(SetParam("maxbaseline", "none"))
	items.append(SetParam("startat", "0"))
	items.append(SetParam("endat", "0"))
	items.append(SetParam("useflag", "true"))
	items.append(SetParam("usemeta", "true"))
	items.append(SetParam("templatebase", "rts_"))

	items.append(SetParam("calcadence", "64"))
	items.append(SetParam("imgcadence", "96"))
	items.append(SetParam("catextent", "20"))
	items.append(SetParam("rtsbin", "rts_gpu"))
	items.append(SetParam("imgcat", "regenerate"))
	items.append(SetParam("calcat", "regenerate"))
	items.append(SetParam("srccat", SRC_CAT))
	items.append(SetParam("calsrcs", "25"))
	items.append(SetParam("imgsrcs", "300"))
	items.append(SetParam("regrid", "false"))
	items.append(SetParam("regridmethod", "3"))
	items.append(SetParam("makepsf", "false"))
	items.append(SetParam("makestokes", "true"))
	items.append(SetParam("removeinstpols", "true"))
	
	logger.info("Reading input file: %s" %(fname))
	obsids = glob.glob("1?????????")
	logger.info("Found %d snapshots in current directory" %(len(obsids)))
	for line in open(fname, "rt"):
		if line[0] == "#":
			continue
		if ":" in line:
			data = line.split(":")
		else:
			data = [line]
		if len(data) < 2:
			if data[0].lower().strip() == "stop":
				break
			logger.warning("Parameters missing: %s" %(line))
			continue
		itemType = data[0].lower().strip()
		pdata = data[1].strip().split()
		if len(pdata) < 1:
			logger.warning("Parameters missing: %s" %(line))
			continue
		if itemType in ["image", "accumulate", "calibrate", "clean", "uvdump", "reflag"]:
			# Check if a range of obsids is specified
			ids = []
			calIds = []
			calSources = []
			if itemType in ["image", "accumulate", "uvdump"]:
				# Check if an ra/dec is supplied to force imaging centre
				targetArgs = extractArgs(pdata[0])
				if len(targetArgs) > 2:
					rahrs = targetArgs[1]
					decdeg = targetArgs[2]
				else:
					rahrs = None
					decdeg = None
				# keep only the id information
				pdata[0] = targetArgs[0]
			rdata = pdata[0].split("-")
			# Add just a single obsid
			if len(rdata) == 1:
				ids.append(pdata[0])
			else:
				# Add a range of snapshots if the obsids are within the specified range
				minid = int(rdata[0].strip())
				maxid = int(rdata[1].strip())
				for obsid in obsids:
					if int(obsid) < minid:
						continue
					if int(obsid) > maxid:
						continue
					ids.append(obsid)
			ids.sort()
			for taskid in ids:
				if itemType == "image":
					if len(pdata) != 3:
						logger.warning("Confused: %s" %(line))
						continue
					items.append(Image(taskid, pdata[1], pdata[2], rahrs, decdeg))
				elif itemType == "clean":
					items.append(Clean(taskid))
				elif itemType == "uvdump":
					if len(pdata) != 3:
						logger.warning("Confused: %s" %(line))
						continue
					items.append(UVdump(taskid, pdata[1], pdata[2], rahrs, decdeg))
				elif itemType == "accumulate":
					if len(pdata) != 3:
						logger.warning("Confused: %s" %(line))
						continue
					items.append(Accumulate(taskid, pdata[1], pdata[2], rahrs, decdeg))
				elif itemType == "calibrate":
					items.append(Calibrate(taskid))
				elif itemType == "fetchmetadata":
					items.append(FetchMetadata(taskid))
#				elif itemType == "flag":
#					items.append(Flag(taskid, " ".join(pdata[1:])))
				elif itemType == "reflag":
					items.append(Reflag(taskid, pdata[1]))
				elif itemType == "pairreflag":
					items.append(PairReflag(taskid, pdata[1], pdata[2]))
				else:
					logger.warning("Confused: %s" %(line))
		elif itemType == "fetchdata":
			items.append(FetchData(pdata[0]))
			# Need to add this to the list of obsids since additional processing may be applied to this obsid.
			obsids.append(pdata[0])
		else:
			items.append(SetParam(data[0].lower().strip(), data[1].lower().strip()))
	return items, parameters

# Set up some logging
if len(sys.argv) < 2:
        sys.exit("Missing input file parameter")

data = sys.argv[1].split(".")
logfile = "%s.log" %(".".join(data[:-1]))
logger = logging.getLogger("rtsprocess")
hdlr = logging.FileHandler(logfile)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

items, params = readProcessFile(sys.argv[1])
for i in items:
	i.execute(params)
