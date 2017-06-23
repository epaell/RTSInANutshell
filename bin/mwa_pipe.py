#!/usr/bin/env python

"""
mwa_pipe.py

python classes/routines to implement simple routines for running the RTS

$Rev$:     Revision of last commit
$Author$:  Author of last commit
$Date$:    Date of last commit


"""
import errno
import os
import sys
import glob
import logging
import time
import numpy as np
from astropy.io import fits
import datetime
import pytz
from dateutil import parser
import json

try:
	# For Python 3.0 and later
	from urllib.request import urlopen
	from urllib.parse import urlencode
	from urllib.error import HTTPError
	from urllib.error import URLError
except ImportError:
	# Fall back to Python 2's urllib2
	from urllib import urlencode
	from urllib2 import urlopen
	from urllib2 import HTTPError
	from urllib2 import URLError

	  
def read_metafits(path):
	"""read header parameters from the metafits file.

	Args:
		path : path in which metafits may be found.

	Returns:
		a dictionary containing all of the header information referenced by header keyword.
	
	Raises:
	
	"""
	header = {}
	radeg = 0.0
	dec_deg = 0.0
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

def get_mwaf_chans(f):
	"""read the number of channels in an mwaf flag file.

	Args:
		path : the mwaf file to read.

	Returns:
		the number of channels in the mwaf file.
	
	Raises:
	
	"""
	hdu = fits.open(f)
	nchan = int(hdu[0].header['NCHANS'])
	hdu.close()
	return nchan

def read_mwaf_data(obs_id):
	"""read the flagging information for all coarse channels.

	Args:
		obs_id : the obs_id for which flagging data will be read.

	Returns:
		an array (nCoarse x nchan) containing the number of visibilities flagged in each channel of each coarse channel.
		the number of channels per coarse channel.

	Raises:

	"""
	COARSE_CHANNELS = 24
	# Read flagging information for obs_id
	mwaf_files = glob.glob("%s_*.mwaf" %(obs_id))
	assert (len(mwaf_files) == COARSE_CHANNELS), "Couldn't find 24 MWAF files."
	nchan = get_mwaf_chans(mwaf_files[0])
	nflagged = np.zeros((COARSE_CHANNELS, nchan), dtype = np.int)

	for cc in range(COARSE_CHANNELS):
		in_file = "%s_%02d.mwaf" %(obs_id, cc + 1)
		hdu = fits.open(in_file)

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

def update_mwaf_data(obs_id, flagged, threshold):
	"""update the flagging information for all coarse channels. Any channels with the number of flagged
		visibilities exceeding the median by a specified threshold will be flagged entirely.

	Args:
		obs_id : the obs_id for which flagging data will be updated.
		flagged : the flagged channels read using read_mwaf_data.
		threshold : the factor above the median number of flagged visibilities at which the entire channel is flagged.

	Returns:
		number of flagged channels
		total number of channels
	
	Raises:

	"""
	COARSE_CHANNELS = 24
	mwaf_files = glob.glob("%s_*.mwaf" %(obs_id))
	if(len(mwaf_files) != COARSE_CHANNELS):
		sys.exit(-1)
	nchan = get_mwaf_chans(mwaf_files[0])

	base = np.median(flagged)
	bad = np.where(flagged > threshold*base)
#	print("Flagging %d channels out of %d" %(len(bad[1]), nchan * COARSE_CHANNELS))
	nbad = len(bad[1])
	for cc in range(COARSE_CHANNELS):
		in_file = "%s_%02d.mwaf" %(obs_id, cc + 1)
#		logger.info("Updating %s:" %(in_file))
		hdu = fits.open(in_file, mode="update")

		n_reflag = 0
		for i in range(nchan):
			if flagged[cc][i] > base * threshold:
#				logger.info(i)
				hdu[1].header['REFLG_%02d' % n_reflag] = i
				n_reflag += 1

		hdu.close()
	return len(bad[1]), nchan * COARSE_CHANNELS

def remove_pols(pols = ["XX", "YY", "XYim", "XYre"]):
	for pol in pols:
		os.system("find -name \"*%s.fits\" -type f -print0 | xargs -0 munlink" %(pol))
	# Try to get a list of all of the instrumental XX images
#	flist = glob.glob("2*_*%s.fits" %(pols[0]))
#	while len(flist) > 0:
#		# Get the frequency of the first item in the list
#		freq0 = flist[0].split("_")[2]
#		# Find all time stamps available for that frequency (need to minimize the number of files being removed at any time)
#		flist = glob.glob("2*%s*%s.fits" %(freq0, pols[0]))
#		for fname in flist:
#			tindex = fname.split("_")[-1]
#			os.system("rm 2*%s" %(tindex))
#			os.system("rm 2*%s" %(tindex.replace(pols[0], pols[1])))
#			os.system("rm 2*%s" %(tindex.replace(pols[0], pols[2])))
#			os.system("rm 2*%s" %(tindex.replace(pols[0], pols[3])))
	# Remove any files
#	os.system("rm 2*%s.fits" %(pols[1]))
#	os.system("rm 2*%s.fits" %(pols[2]))
#	os.system("rm 2*%s.fits" %(pols[3]))
	return

def move_images(pols, dest):
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

class MWAPipe:
	def __init__(self, pipe_name):
		# Set up the logger
		self.CAT_PATH = "/home/elenc/catalog/"
		self.SRC_CAT = "rtscat.txt"
		log_file = "%s.log" %(pipe_name)
		self.logger = logging.getLogger("mwapipe")
		hdlr = logging.FileHandler(log_file)
		formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
		hdlr.setFormatter(formatter)
		self.logger.addHandler(hdlr)
		self.logger.setLevel(logging.INFO)
		
		cat_path = os.environ.get('RTS_CAT_PATH')
		if cat_path == None:
			self.logger.warning("RTS_CAT_PATH not set, using %s as default" %(self.CAT_PATH))
			cat_path = self.CAT_PATH
		self.cat_path = cat_path
		self.data_path = None			# Path where data is stored
		self.weighting = "natural"
		self.robustness = 2.0
		self.do_accumulate = False
		self.work_dir = os.getcwd()
		self.array = "128T"
		self.field_size = 2.0
		self.nside = 2048
		self.max_subbands = -1
		self.fscrunch = 4

		self.adjust_cal_params = True
		self.ncalibrators = 1
		self.npeel = 0
		self.niono_cal = 0
		self.update_cal = False
		self.short_taper = None
		self.long_taper = None
		self.min_baseline = None
		self.max_baseline = None
		self.start_at = 0
		self.end_at = 0
		self.use_flag = True
		self.use_meta = True
		self.template_base = "rts_"

		self.cal_cadence = 64
		self.img_cadence = 96
		self.cat_extent = 20.0
		self.rts_bin = "rts_gpu"
		self.img_cat = None
		self.cal_cat = None
		self.src_cat = self.SRC_CAT
		self.cal_srcs = 25
		self.img_srcs = 300
		self.regrid = False
		self.regrid_method = 3
		self.make_psf = False
		self.make_stokes = True
		self.remove_inst_pols = True

	def get_data_path(self, obs_id = None):
		if self.data_path == None:
			if obs_id == None:
				return self.work_dir
			else:
				return "%s/%s" %(self.work_dir, obs_id)
		if obs_id == None:
			return self.data_path
		return "%s/%s" %(self.data_path, obs_id)

	def move_to_target(self, obs_id = None):
		if obs_id == None:
			work_path = self.work_dir
		else:
			work_path = "%s/%s" %(self.work_dir, obs_id)
		# Make sure the destination exists
		os.system("mkdir -p %s" %(work_path))
		# Change directory
		os.chdir(work_path)
		
	def move_to_data(self, obs_id = None):
		data_path = self.get_data_path(obs_id)
		assert (os.path.exists(data_path) == True), "Data path does not exist: %s" %(data_path)
		os.chdir(data_path)
		
	def return_home(self):
		# Change directory
		os.chdir(self.work_dir)
		
	def _generate_RTS_input_files(self, obs_id, basename, template, ra_hrs, dec_deg, mode = None):
		import ephem
		from mwapy import ephem_utils

		if not (self.array == '32T' or self.array == '128T'):
			self.logger.error("Array Parameter should be either '32T' or '128T'")
			return
	
		data_path = self.get_data_path(obs_id)
		obs_path = "%s/%s" %(self.work_dir, obs_id)
		
		# Set the time.
		self.logger.info("obs_id=%s" %(obs_id))
		obs_time = ephem_utils.MWATime(gpstime = float(obs_id))
		lst_hrs = float(obs_time.LST) / 15.0
		jd = obs_time.MJD + 2400000.5
		self.logger.info("lst=%.4f hr" %(lst_hrs))
		self.logger.info("jd=%.5f" %(jd))
		self.logger.info("Reading metafits information for %s" %(data_path))
		header = read_metafits(data_path)
		bandwidth = header["BANDWIDTH"] / 24
		if header["HA_HRS"] > 12.0:
			header["HA_HRS"] -= 24.0

		# How many GPU boxes are there?
		gpu_files = glob.glob(data_path + '/*gpubox*00.fits')  

		if(len(gpu_files) > 0):
			band_list = [int(filename[-10:-8]) for filename in gpu_files]
		else:
			uv_file_list = glob.glob(data_path + '/*uvfits')
			band_list = [int(filename[-9:-7]) for filename in uv_file_list]   
		band_list.sort()
		if self.max_subbands != -1:
			band_list = band_list[:self.max_subbands]	# Limit the number of subbands to MaxSubBands
		subband_string = ''
		for band in band_list:
			subband_string = subband_string + str(band) + ','
		subband_string = subband_string[:-1]

		corr_dump_time = header["INT_TIME"]

		# Generate a channel flag file to remove centre and edge channels
		template_file = open(template)
		set_subbands = 0

		# Write a new .in file
		outName = "rts_%s.in" %(basename)
		outFile = open(outName, "w+")
	
		corr_dumps_per_cadence = -1.0
		integration_bins = 0.0
		if mode == "calibrate":
			outFile.write("generateDIjones=1\n")
			outFile.write("useStoredCalibrationFiles=0\n")
			imaging_cadence = self.cal_cadence
			integration_bins = 6 - int(np.log2(corr_dump_time))
		else:
			imaging_cadence = self.img_cadence
		
			integration_bins = 4 - int(np.log2(corr_dump_time))
		corr_dumps_per_cadence = (int)(imaging_cadence / corr_dump_time)
		if corr_dumps_per_cadence == 1:
			integration_bins = 1
		elif corr_dumps_per_cadence == 2:
			integration_bins = 2
		elif corr_dumps_per_cadence < 8:
			integration_bins = 3
		if mode == "uv":
			outFile.write("generateDIjones=0\n")
			outFile.write("useStoredCalibrationFiles=1\n")
		if mode == "image":
			outFile.write("generateDIjones=0\n")
			outFile.write("useStoredCalibrationFiles=1\n")
			if self.weighting == "uniform":
				outFile.write("imgUniformWeighting=1\n")
			elif self.weighting == "robust":
				outFile.write("imgRobustWeighting=1\n")
				outFile.write("imgRobustnessParameter=%f\n" %(self.robustness))
			if self.do_accumulate == True:
				outFile.write("UseAccumUVsampling=1\n")
			if self.regrid == True:
				outFile.write("MakeImage=\n")
				outFile.write("FieldOfViewDegrees=%f\n" %(self.field_size))
				outFile.write("RegridMethod=%d\n" %(self.regrid_method))
				outFile.write("DoRegriddingWithProjection=%d\n" %(self.nside))
			else:
				outFile.write("FieldOfViewDegrees=%f\n" %(self.field_size))
				if self.make_stokes == True:
					outFile.write("MakeStokesSnapshots=\n")
				else:
					outFile.write("MakeImage=\n")
			if self.make_psf == True:
				outFile.write("ImagePSF=1\n")
		if mode == "accumulate":
			outFile.write("generateDIjones=0\n")
			outFile.write("useStoredCalibrationFiles=1\n")
			outFile.write("DoUVsamplingAccum=1\n")
			outFile.write("FieldOfViewDegrees=%f\n" %(self.field_size))
			outFile.write("MakeStokesSnapshots=\n")
		if mode == "image" or mode == "accumulate":
			if self.short_taper != None:
				outFile.write("imgShortBaselineTaper=%f\n" %(self.short_taper))
			if self.long_taper != None:
				outFile.write("imgLongBaselineTaper=%f\n" %(self.long_taper))
			if self.min_baseline != None:
				outFile.write("imgBaselineMin=%f\n" %(self.min_baseline))
			if self.max_baseline != None:
				outFile.write("imgBaselineMax=%f\n" %(self.max_baseline))

		self.logger.info("Imaging cadence = %f" %(imaging_cadence))
		scan_time = header["N_SCANS"] * header["INT_TIME"] - header["INT_TIME"] * 4
		self.logger.info("scan_time = %d" %(scan_time))
		n_iterations = int(scan_time / imaging_cadence)
		self.logger.info("n_iterations = %d" %(n_iterations))

		template_file = open(template)
		for line in template_file:
			new_line = line
			if (self.start_at != 0) and (line.find('StartProcessingAt') == 0):
				new_line = "StartProcessingAt=%d\n" %(self.start_at)
			if (self.start_at != 0) and (line.find('StartIntegrationAt') == 0):
				new_line = "StartIntegrationAt=%d\n" %(self.start_at)
			if line.startswith('BaseFilename'):
				new_line = "BaseFilename=%s/*_gpubox\n" %(data_path)
				if self.use_flag == True:
					new_line += "ImportCotterFlags=1\nImportCotterBasename=%s/%s\n\n" %(obs_path, obs_id)
				if self.use_meta == True:
					new_line += "ReadMetafitsFile=1\nMetafitsFilename=%s/%s\n" %(data_path, obs_id)
			if line.startswith('CorrDumpTime'):
				new_line = line.replace(line[len('CorrDumpTime='):],"%.1f" %(header["INT_TIME"]) + '\n')
			if line.startswith('FscrunchChan'):
				new_line = "FscrunchChan=%d\n" %(self.fscrunch)
			if line.startswith('ChannelBandwidth'):
				new_line = line.replace(line[len('ChannelBandwidth='):],"%.3f" %(header["BANDWIDTH"]/header["N_CHANS"]) + '\n')
			if line.startswith('NumberOfChannels'):
				new_line = line.replace(line[len('NumberOfChannels='):],str(header["N_CHANS"]/24) + '\n')
			if line.startswith('ObservationFrequencyBase'):
				new_line = line.replace(line[len('ObservationFrequencyBase='):],str(header["FREQCENT"] - header["BANDWIDTH"] / 2.0 - 0.02) + '\n')
			if line.startswith('ObservationPointCentreHA'):
				new_line = line.replace(line[len('ObservationPointCentreHA='):],str(header["HA_HRS"]) + '\n')
			if line.startswith('ObservationPointCentreDec'):
				new_line = line.replace(line[len('ObservationPointCentreDec='):],str(header["DEC_DEGS"]) + '\n')
			if 'SourceCatalogueFile' in line:
				if mode == "calibrate":
					self.logger.info("Setting catalogue file for self calibration")
					new_line = "SourceCatalogueFile=catalogue.txt\n"
				else:
					self.logger.info("Setting catalogue file for imaging")
					new_line = "SourceCatalogueFile=catalogue_imaging.txt\n"

			if self.adjust_cal_params == True:
				if 'NumberOfSourcesToPeel' in line:
						new_line = "NumberOfSourcesToPeel=%s\n" %(self.npeel)
				if 'NumberOfIonoCalibrators' in line:
						new_line = "NumberOfIonoCalibrators=%s\n" %(self.niono_cal)
				if 'UpdateCalibratorAmplitudes' in line:
						if self.update_cal == True:
							new_line = "UpdateCalibratorAmplitudes=1\n"
						else:
							new_line = "UpdateCalibratorAmplitudes=0\n"
			if line.startswith('ObservationImageCentreRA'):
				if ra_hrs == None:
					new_line = line.replace(line[len('ObservationImageCentreRA='):],str(header["RA_HRS"]) + '\n')
				else:
					if ra_hrs == "zenith":
						new_line = "ObservationImageCentreRA=%f\n" %(lst_hrs)
					else:
						new_line = "ObservationImageCentreRA=%f\n" %(ra_hrs)
			if line.startswith('ObservationImageCentreDec'):
				if dec_deg == None:
					new_line = line.replace(line[len('ObservationImageCentreDec='):],str(header["DEC_DEGS"]) + '\n')
				else:
					if dec_deg == "zenith":
						new_line = "ObservationImageCentreDec=%f\n" %(-26.0 - 42.0/60.0 - 11.95/3600.0)
					else:
						new_line = "ObservationImageCentreDec=%f\n" %(dec_deg)
			if line.startswith('CorrDumpsPerCadence'):
				new_line = "CorrDumpsPerCadence=%d\n" %(corr_dumps_per_cadence)
			if line.startswith('NumberOfIntegrationBins'):
				new_line = "NumberOfIntegrationBins=%d\n" %(integration_bins)
			if line.startswith('NumberOfIterations'):
				if self.end_at != 0:
					new_line = "NumberOfIterations=%d\n" %(self.end_at)
				else:
					if(self.array=='32T'):
						# corr_dumps_per_cadence is currently hard set to 4
						new_line = line.replace(line[len('NumberOfIterations='):],str(header["N_SCANS"] / 4) + '\n')
					else:
						new_line = line.replace(line[len('NumberOfIterations='):],str(n_iterations) + '\n')
			if line.startswith('SubBandIDs='):
				new_line = line.replace(line[len('SubBandIDs='):],subband_string)
				set_subbands = 1
			if line.startswith('ObservationTimeBase='):
				new_line = "ObservationTimeBase=%.5f\n" %(jd)
	#			new_line = line.replace(line[len('ObservationTimeBase='):], "%.5f" %(jd))	
			outFile.write(new_line)
		if(self.array=='128T'):
			if(set_subbands == 0):
				outFile.write('SubBandIDs='+subband_string+'\n')
		outFile.close()
		template_file.close()
		return len(band_list)

	def clean(self, obs_id):
		"""removes all temporary files from the specified target obs_id.

		Args:
			obs_id : the obs_id where all temporary files will be removed.

		Returns:
	
		Raises:
	
		"""
		self.move_to_target(obs_id)
		# Remove old files
		self.logger.info("Cleaning up %s/%s" %(self.work_dir, obs_id))
		# Remove instrumental polarisation images
		remove_pols()
		# Remove Stokes images
		remove_pols(["I", "Q", "U", "V"])
		# Remove any other miscellaneous files
		os.system("rm -fr *.log 2*.fits *.btr uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits *.dat cat*.txt rts.in")
		# Return to working directory
		self.return_home()

	def reflag(self, obs_id, threshold):
		"""read flagging information and completely flag channels that already have excessive flagging.

		Args:
			obs_id : the obs_id where excessively flagged channels will be flagged entirely.
			threshold : the factor above the median number of flagged visibilities at which the entire channel is flagged.

		Returns:
	
		Raises:
	
		"""
		self.logger.info("Reflagging %s : threshold=%f" %(obs_id, threshold))
		# Work out where the data resides
		data_path = self.get_data_path(obs_id)
		# Move to the work directory
		self.move_to_target(obs_id)
		# Remove any old mwaf files that may exist in the work path
		os.system("rm -fr *.mwaf")
		# Extract mwaf files from the zip file downloaded from the archive
		zip_file = "%s/%s_flags.zip" %(data_path, obs_id)
		if os.path.isfile(zip_file) == False:
			self.logger.warning("Cotter flag zip file for obs_id %s does not exist. WARNING: not using Cotter flags." %(obs_id))
		os.system("unzip -o -j %s" %(zip_file))

		# Read existing flagging information
		try:
			flagged, nchan = read_mwaf_data(obs_id)
		except AssertionError as e:
			self.logger.error(e.args[0])
			return
		# Flag excessively flagged channels
		nbad, ntotal = update_mwaf_data(obs_id, flagged, threshold)
		self.logger.info("Flagging %d channels out of %d" %(nbad, ntotal))
		self.logger.info("Finished reflagging %s" %(obs_id))
		# Return to working directory
		self.return_home()

	def pair_reflag(self, obs_id1, obs_id2, threshold):
		"""read flagging information for a pair of target obs_ids and completely flag channels that already
			have excessive flagging. Flagging in both targets will be made consistent so that the same channels
			are flagged in each obs_id.

		Args:
			obs_id1 : the first obs_id where excessively flagged channels will be flagged entirely.
			obs_id2 : the second obs_id where excessively flagged channels will be flagged entirely.
			threshold : the factor above the median number of flagged visibilities at which the entire channel is flagged.

		Returns:
	
		Raises:
	
		"""
		self.logger.info("Reflagging Pair %s+%s : threshold=%f" %(obs_id1, obs_id2, threshold))
		# Work out where the data resides
		data_path = self.get_data_path(obs_id1)
		# Move to the work directory
		self.move_to_target(obs_id1)
		# Remove any old mwaf files that may exist in the work path
		os.system("rm -fr *.mwaf")
		# Extract mwaf files from the zip file downloaded from the archive
		zip_file = "%s/%s_flags.zip" %(data_path, obs_id1)
		if os.path.isfile(zip_file) == False:
			self.logger.warning("Cotter flag zip file for obs_id %s does not exist. WARNING: not using Cotter flags." %(obs_id1))
		os.system("unzip -o -j %s" %(zip_file))

		# Read existing flagging information
		try:
			flagged1, nchan1 = read_mwaf_data(obs_id1)
		except AssertionError as e:
			self.logger.error(e.args[0])
			return

		# Work out where the data for the second obs_id resides
		data_path = self.get_data_path(obs_id2)
		# Move to the work directory for that obs_id
		self.move_to_target(obs_id2)
		# Remove any old mwaf files that may exist in the work path
		os.system("rm -fr *.mwaf")
		# Extract mwaf files from the zip file downloaded from the archive
		zip_file = "%s/%s_flags.zip" %(data_path, obs_id2)
		if os.path.isfile(zip_file) == False:
			self.logger.warning("Cotter flag zip file for obs_id %s does not exist. WARNING: not using Cotter flags." %(obs_id2))
		os.system("unzip -o -j %s" %(zip_file))

		# Read existing flagging information
		try:
			flagged2, nchan2 = read_mwaf_data(obs_id2)
		except AssertionError as e:
			self.logger.error(e.args[0])
			return
		if nchan1 != nchan2:
			self.logger.error("Channels in specified obs_ids do not match %s=%d; %s=%d" %(obs_id1, nchan1, obs_id2, nchan2))
			return
		flagged = flagged1 + flagged2
		nchan = nchan1
		self.move_to_target(obs_id1)
		# Flag excessively flagged channels in the first target
		nbad, ntotal = update_mwaf_data(obs_id1, flagged, threshold)
		self.logger.info("Flagging %d channels out of %d" %(nbad, ntotal))
		self.move_to_target(obs_id2)
		# Flag excessively flagged channels in the second target
		nbad, ntotal = update_mwaf_data(obs_id2, flagged, threshold)
		self.logger.info("Flagging %d channels out of %d" %(nbad, ntotal))
		self.logger.info("Finished reflagging %s+%s" %(obs_id1, obs_id2))
		# Return to working directory
		self.return_home()

	def fetch_data(self, obs_id):
		"""fetch visilibity, flagging and metadata files from the archive for the specified obs_id.

		Args:
			obs_id : the obs_id for which data will be retrieved from the archive.

		Returns:
	
		Raises:
	
		"""
		self.logger.info("Fetch data for %s" %(obs_id))
		self.move_to_data()
		os.system("obsdownload.py -o %s" %(obs_id))
		self.logger.info("Finished fetching data for %s" %(obs_id))
		self.return_home()

	def fetch_metadata(self, obs_id):
		"""regenerate metadata files (rather than fetching from the archive) for the specified obs_id.

		Args:
			obs_id : the obs_id for which metadata will be generated.

		Returns:
	
		Raises:
	
		"""
		self.logger.info("Fetch metadata for %s" %(obs_id))
		# Move to the target directory
		self.move_to_data(obs_id)
		os.system("wget -O %s_metafits_ppds.fits http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=%s" %(obs_id, obs_id))
		# Return to working directory
		self.return_home()
	
	def calibrate(self, obs_id):
		self.logger.info("Calibrate %s" %(obs_id))
		self.move_to_target(obs_id)

		# Remove old files
		self.logger.info("Removing temporary files in %s/%s" %(self.work_dir, obs_id))
		os.system("rm -fr *.log 2*.fits *.btr uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits *.dat")

		# Use selfcal of field to do calibration
		self.logger.info("Generating initial calibration with %s" %(obs_id))
		if self.cal_cat == None:
			metafits_path = "%s/%s_metafits_ppds.fits" %(self.get_data_path(obs_id), obs_id)
			if os.path.exists(metafits_path) == True:
				self.logger.info("Using %s as input" %(metafits_path))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -m %s -c %f -s %s/%s -n %s" %(metafits_path, self.cat_extent, self.cat_path, self.src_cat, self.cal_srcs))
			else:
				self.logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(obs_id))
			os.system("mv *patch*.txt catalogue.txt")
		else:
			# Re-use a catalogue generated for the specified obs_id
			os.system("cp %s/%s/catalogue.txt ." %(self.work_dir, self.cal_cat))

		# Generate an RTS input file based on the template
		self.logger.info("Generating input files %s" %(obs_id))
		ncoarse = self._generate_RTS_input_files(obs_id, "cal", "%s/%scal.in" %(self.work_dir, self.template_base), None, None, "calibrate")
		# Copy over the default flag files
#				os.system("cp %s/fl*.txt ." %(self.work_dir))

		# Do the calibration
		self.logger.info("Performing calibration for %s" %(obs_id))
		os.system("aprun -N 1 -n %d %s rts_cal.in" %(ncoarse + 1, self.rts_bin))
			
		self.logger.info("Finished calibrating %s" %(obs_id))
		# Return to working directory
		self.return_home()

	def image(self, obs_id, dest_image_path, cal_id = None, ra_hrs = None, dec_deg = None):
		self.logger.info("Imaging %s" %(obs_id))
		self.logger.info("Using calibration from %s" %(cal_id))
		if self.weighting == "robust":
			self.logger.info("Weighting = %s (robustness=%.1f)" %(self.weighting, self.robustness))
		else:
			self.logger.info("Weighting = %s" %(self.weighting))

		if ra_hrs != None and dec_deg != None:
			self.logger.info("Using (%f, %f) as imaging centre" %(ra_hrs, dec_deg))
		# Move to the target directory
		self.logger.info("Cleaning up %s/%s" %(self.work_dir, obs_id))
		self.move_to_target(obs_id)
		# Remove old files
		remove_pols()
		remove_pols(["I", "Q", "U", "V"])
		os.system("rm -fr *.log *.btr pee*.txt rest*.txt cor*.fits int*.fits")
		# Move to final cal to see if there is anything there
		if cal_id == None:
			self.move_to_target(obs_id)
			# Use calibration already in target
		else:
			# Use calibration from specified obs_id
			self.move_to_target(cal_id)
			if cal_id != obs_id:
				# Copy the calibration solutions to the target
				self.logger.info("Copying calibration files from %s to %s" %(cal_id, obs_id))
				# copy across the calibration data ... unless we are imaging the calibrator.
				os.system("cp *.dat %s/%s" %(self.work_dir, obs_id))

		# Check if a calibration exists
		if len(glob.glob("*.dat")) == 0:
			self.logger.warning("No calibration data found!")
			return

		self.move_to_target(obs_id)
		# Generate the catalogue for imaging.
		if self.img_cat == None:
			metafits_path = "%s/%s_metafits_ppds.fits" %(self.get_data_path(obs_id), obs_id)
			if os.path.exists(metafits_path) == True:
				self.logger.info("Using %s as input" %(metafits_path))
				self.logger.info("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s -c %f -s %s/%s -n %s" %(metafits_path, self.cat_extent, self.cat_path, self.src_cat, self.img_srcs))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s -c %f -s %s/%s -n %s" %(metafits_path, self.cat_extent, self.cat_path, self.src_cat, self.img_srcs))
			else:
				self.logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(obs_id))
			os.system("mv *peel*.txt catalogue_imaging.txt")
		else:
			self.logger.info("Using catalogue from %s" %(self.img_cat))
			os.system("cp %s/%s/catalogue_imaging.txt ." %(self.work_dir, self.img_cat))
			
		# Generate an RTS input file based on the template
		self.logger.info("Generating RTS imaging input files for %s" %(obs_id))
		ncoarse = self._generate_RTS_input_files(obs_id, "img", "%s/%simg.in" %(self.work_dir, self.template_base), ra_hrs, dec_deg, "image")
		# Do the imaging
		self.logger.info("Imaging %s" %(obs_id))
		os.system("aprun -N 1 -n %d %s rts_img.in" %(ncoarse + 1, self.rts_bin))
		if self.remove_inst_pols == True:
			remove_pols()
		self.logger.info("Moving image files to destination path: %s" %(dest_image_path))
		os.system("mkdir -p %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		move_images(["XX", "YY", "XYim", "XYre"], "%s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		move_images(["I", "Q", "U", "V"], "%s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("mv peeled_sources_*.txt %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("mv restore_*.txt %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("mv rts_*.log %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("cp *.dat %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("cp *.in %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		self.logger.info("Finished imaging %s" %(obs_id))
		# Return to working directory
		self.return_home()

	def accumulate(self, obs_id, dest_image_path, cal_id = None, ra_hrs = None, dec_deg = None):
		self.logger.info("Accumulating weights for %s" %(obs_id))
		self.logger.info("Using calibration from %s" %(cal_id))

		if ra_hrs != None and dec_deg != None:
			self.logger.info("Using (%f, %f) as imaging centre" %(ra_hrs, dec_deg))
		# Move to the target directory
		self.move_to_target(obs_id)
		self.logger.info("Cleaning up %s/%s" %(self.work_dir, obs_id))
		# Remove old files
		os.system("rm -fr *.log 2*.fits *.btr uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits")
		# Move to final cal to see if there is anything there
		if cal_id == None:
			self.move_to_target(obs_id)
			# Use calibration already in target
		else:
			self.move_to_target(cal_id)
			# Use calibration from specified obs_id
			if cal_id != obs_id:
				# Copy the calibration solutions to the target
				self.logger.info("Copying calibration files from %s to %s" %(cal_id, obs_id))
				# copy across the calibration data ... unless we are imaging the calibrator.
				os.system("cp *.dat %s/%s" %(self.work_dir, obs_id))

		# Check if a calibration exists
		if len(glob.glob("*.dat")) == 0:
			self.logger.warning("No calibration data found!")
			return

		self.move_to_target(obs_id)
		# Generate the catalogue for imaging.
		if self.img_cat == None:
			metafits_path = "%s/%s_metafits_ppds.fits" %(self.get_data_path(obs_id), obs_id)
			if os.path.exists(metafits_path) == True:
				self.logger.info("Using %s as input" %(metafits_path))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s -c %f -s %s/%s -n %s" %(metafits_path, self.cat_extent, self.cat_path, self.src_cat, self.img_srcs))
			else:
				self.logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(obs_id))
			os.system("mv *peel*.txt catalogue_imaging.txt")
		else:
			self.logger.info("Using catalogue from %s" %(self.img_cat))
			os.system("cp %s/%s/catalogue_imaging.txt ." %(self.work_dir, self.img_cat))
		
		# Generate an RTS input file based on the template
		self.logger.info("Generating RTS input files to accumulate weights for %s" %(obs_id))
		ncoarse = self._generate_RTS_input_files(obs_id, "acc", "%s/%simg.in" %(self.work_dir, self.template_base), ra_hrs, dec_deg, "accumulate")
		# Generate the weights file
		self.logger.info("Accumulate weights for %s" %(obs_id))
		os.system("aprun -N 1 -n %d %s rts_acc.in" %(ncoarse + 1, self.rts_bin))

		os.system("mv rts_*.log %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("cp *.dat %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		self.logger.info("Finished accumulating weights for %s" %(obs_id))
		# Return to working directory
		self.return_home()

	def uv_dump(self, obs_id, dest_image_path, cal_id = None, ra_hrs = None, dec_deg = None):
		self.logger.info("UV Dump for %s" %(obs_id))
		self.logger.info("Using calibration from %s" %(cal_id))

		if ra_hrs != None and dec_deg != None:
			self.logger.info("Using (%f, %f) as imaging centre" %(ra_hrs, dec_deg))
		# Move to the target directory
		self.move_to_target(obs_id)
		self.logger.info("Cleaning up %s/%s" %(self.work_dir, obs_id))
		# Remove old files
		os.system("rm -fr *.log 2*.fits *.btr uvdump*.fits uvbeam*.fits pee*.txt rest*.txt cor*.fits int*.fits")
		# Move to final cal to see if there is anything there
		if cal_id == None:
			self.move_to_target(obs_id)
			# Use calibration already in target
		else:
			# Use calibration from specified obs_id
			self.move_to_target(cal_id)
			if cal_id != obs_id:
				# Copy the calibration solutions to the target
				self.logger.info("Copying calibration files from %s to %s" %(cal_id, obs_id))
				# copy across the calibration data ... unless we are imaging the calibrator.
				os.system("cp *.dat %s/%s" %(self.work_dir, obs_id))

		# Check if a calibration exists
		if len(glob.glob("*.dat")) == 0:
			self.logger.warning("No calibration data found!")
			return

		self.move_to_target(obs_id)
		# Generate the catalogue for imaging.
		if self.img_cat == None:
			metafits_path = "%s/%s_metafits_ppds.fits" %(self.get_data_path(obs_id), obs_id)
			if os.path.exists(metafits_path) == True:
				self.logger.info("Using %s as input" %(metafits_path))
				os.system("python /group/mwaops/CODE/bin/srclist_by_beam.py -x -m %s -c %f -s %s/%s -n %s" %(metafits_path, self.cat_extent, self.cat_path, self.src_cat, self.img_srcs))
			else:
				self.logger.error("Unable to find %s_metafits_ppds.fits to generate catalogue against!" %(obs_id))
			os.system("mv *peel*.txt catalogue_imaging.txt")
		else:
			self.logger.info("Using catalogue from %s" %(self.img_cat))
			os.system("cp %s/%s/catalogue_imaging.txt ." %(self.work_dir, self.img_cat))
		
		# Generate an RTS input file based on the template
		self.logger.info("Generating RTS input files to accumulate weights for %s" %(obs_id))
		ncoarse = self._generate_RTS_input_files(obs_id, "uv", "%s/%suv.in" %(self.work_dir, self.template_base), ra_hrs, dec_deg, "uv")
		# Generate the weights file
		self.logger.info("UV Dump for %s" %(obs_id))
		os.system("aprun -N 1 -n %d %s rts_uv.in" %(ncoarse + 1, self.rts_bin))

		os.system("mkdir -p %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("mv rts_*.log %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("cp *.dat %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		os.system("mv uvdump*.uvfits %s/%s/%s" %(self.work_dir, dest_image_path, obs_id))
		self.logger.info("Finished dumping uv for %s" %(obs_id))
		# Return to working directory
		self.return_home()

# Function to call a JSON web service and return a dictionary:
def getmeta(service='obs', params=None):
	"""Given a JSON web service ('obs', find, or 'con') and a set of parameters as
		a Python dictionary, return a Python dictionary containing the result.
	"""
	BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'
	if params:
		data = urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
	else:
		data = ''
	# Validate the service name
	if service.strip().lower() in ['obs', 'find', 'con']:
		service = service.strip().lower()
	else:
		print("invalid service name: %s" % service)
		return
	# Get the data
	try:
		if (sys.version_info > (3, 0)):
			# Python 3 works differently ... of course ;-)
			response = urlopen(BASEURL + service + '?' + data)
			str_response = response.readall().decode('utf-8')
			result = json.loads(str_response)
		else:
			result = json.load(urlopen(BASEURL + service + '?' + data))
	except HTTPError as error:
		print("HTTP error from server: code=%d, response:\n %s" % (error.code, error.read()))
		return
	except URLError as error:
		print("URL or network error: %s" % error.reason)
		return
	# Return the result dictionary
	return result

def get_obs_list(params):
	""" Query the archive to get an observation list meeting the required parameters

		mintime, maxtime: Minimum and maximum values of 'starttime' for the observation.
			Times are INCLUSIVE of the endpoints, eg comparisons use 'starttime>=mintime and startttime<=maxtime'.
		minra, maxra: Minimum and maximum values for ra_pointing for the observation,
		 	in degrees. NOTE - this does NOT handle the 0/360 degree wrap, so if you need a range that spans 0/360 degrees, use two queries and join the results.
		mindec, maxdec: Minimum and maximum values for dec_pointing for the observation, in degrees.
		minel, maxel: Minimum and maximum values for elevation_pointing for the observation, in degrees.
		minaz, maxaz: Minimum and maximum values for azimuth_pointing for the observation, in degrees.
		minsunel, maxsunel: Minimum and maximum values for sun_elevation for the observation, in degrees.
		minsunpd, maxsunpd: Minimum and maximum values for sun_pointing_distance for the observation, in degrees.
		projectid: A string to match against the projectid field for the observation.
		mode: A string to match against the observing mode, eg 'HW_LFILES', 'VOLTAGE_START', etc.
		creator: A string to match against the name of the person creating that observation.
		obsname: A string to match against the observation name.
		contigfreq: Pass 1 to match only observations with 24 contiguous frequency channels,
			or 0 to match only observations with gaps in the frequency channels.
		calibration: Pass 1 to match only calibration observations, or 0 to match only observations
			that are NOT calibrations.
		cenchan: Select only observations with this centre frequency channel.
		anychan: Select only observations that include this channel in their list of 24 frequency channels.
		freq_res: Select only observations with the correlator set to this frequency resolution
			(currently, the only valid numbers are 10, 20, and 40 kHz).
		int_time: Select only observations with the correlator set to this integration time
			(currently, the only valid numbers are 0.5, 1.0, 2.0 and 4.0 seconds).
		future: Pass 1 to match only observations that have not yet happened, or 0 to match
			only observations that are in the past.
		gridpoint: An integer to match against the gridpoint_number column in the schedule_metadata table.
		minfiles: Select only observations with at least this many data files recorded.
		limit: By default, the web service returns a maximum of 100 observations.
			Use this parameter to specify a smaller number, or (with care) a larger one.

		Sample Usage:
			obslist = get_obs_list(params={'obsname': 'UVCeti_121', 'mintime': 1133864600, 'maxtime': 1133864960})
	"""
	obslist = getmeta(service="find", params=params)
	sublist =  [x[0] for x in obslist]
	sublist.sort()
	return sublist

def get_local_obs_list(omin = "0000000000", omax = "9999999999"):
	obs_ids = glob.glob("[0-9]"*10)
	obs_ids.sort()
	return [obs for obs in obs_ids if obs >= omin and obs <= omax]

def get_obs_time(obs_id):
	timestr = getmeta(service="tconv", params={'gpssec':obs_id})
	return parser.parse(timestr)

