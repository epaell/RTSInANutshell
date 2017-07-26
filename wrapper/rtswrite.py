#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 15:46:12 2017
RTSwrite
@author: ianbrown
"""

def write_launch(obsid,fscrunch,img_cad,path,weighting,robustness=(-1),field_size='5',email='none',oversample=5,maxfreq=180):
    launch=['#!/bin/bash -l']
    launch.append('#SBATCH --nodes=1')
    launch.append('#SBATCH --ntasks=1')
    launch.append('#SBATCH --ntasks-per-node=1')
    launch.append('#SBATCH --time=0:20:00')
    launch.append('#SBATCH --partition=gpuq')
    launch.append('#SBATCH --account=mwasci')
    launch.append('#SBATCH --export=NONE')
    launch.append('#SBATCH --output="%s.log"'% obsid)
    launch.append('#SBATCH --open-mode=append')
    if email != 'none':
        launch.append('#SBATCH --mail-type=FAIL')
        launch.append('#SBATCH --mail-user=%s'% email)
    launch.append('#')
    launch.append('#')
    launch.append('wrapper.py --obsid=%s --fscrunch=%d --imgcad=%.1f --path=%s --weighting=%s --robustness=%.1f --size=%.1f --email=%s --oversample=%d --maxfreq=%d'% (obsid,fscrunch,img_cad,path,weighting,robustness,field_size,email,oversample,maxfreq))
    
    #write rtslaunch.sh
    file1 = open('%srts.sh'% obsid,'w')
    file2 = open('%s.log' % obsid,'w')
    for i in range (0,len(launch)):
        lnch=('%s\n'% launch[i])
        file1.write(lnch)
    file2.write(launch[len(launch)-1])
    file1.flush()
    file2.flush()
    file1.close
    file2.close
    return

def write_qload(obsid,email):
    '''
    assigns and writes qload.sh
    '''
    qload=['#!/bin/bash -l']
    qload.append('#SBATCH --nodes=1')
    qload.append('#SBATCH --ntasks=1')
    qload.append('#SBATCH --ntasks-per-node=1')
    qload.append('#SBATCH --time=1:00:00')
    qload.append('#SBATCH --partition=copyq')
    qload.append('#SBATCH --clusters=zeus')
    qload.append('#SBATCH --account=mwasci')
    qload.append('#SBATCH --export=NONE')
    
    if email != 'none':
        qload.append('#SBATCH --mail-type=FAIL')
        qload.append('#SBATCH --mail-user=%s'% email)
        
    qload.append('#SBATCH --output="%s.log"'% obsid)
    qload.append('#SBATCH --open-mode=append')
    qload.append('#')
    qload.append('module load mpifileutils')
    qload.append('#')
    qload.append('rm -f load.log')
    qload.append('mpirun -n 1 ./%sload.py'% obsid)
    
    #write qload.sh
    file1 = open('%sqload.sh'% obsid,'w')
    for i in range (0, len(qload)):
        q=('%s\n'% qload[i])
        file1.write(q) 
    file1.flush()
    file1.close()
    return

def write_load(obsid):
    '''
    assigns values and writes load.py 
    '''
    load=['#!/usr/bin/env python']
    load.append('from mwa_pipe import *')
    load.append('import timeit\n')
    load.append('rts = MWAPipe("load")\n')
    #load.append('# Get a list of obsids ')
    #load.append('obsid = %s'% obsid)
    #load.append('for obsid in obsid:')
    load.append('start = timeit.default_timer()')
    load.append('sys.stdout.write("********************************\\n")')
    load.append('sys.stdout.write("Begin fetching observation data \\n")')
    load.append('sys.stdout.write("********************************\\n")')
    load.append('sys.stdout.flush()')
    load.append('# Fetch the data from the archive')
    load.append('rts.fetch_data(%s)'% obsid)
    load.append('# Reflag the observation')
    load.append('rts.reflag(%s, 1.3)'% obsid)
    load.append('sys.stdout.write("***********************************\\n")')
    load.append('sys.stdout.write("Finished fetching observation data \\n")')
    load.append('sys.stdout.write("***********************************\\n")')
    load.append('print("total time to fetch data:%f")%(timeit.default_timer()-start)')    
    
    #write load.py
    file2 = open('%sload.py'% obsid,'w')
    for i in range (0,len(load)):
        l=('%s\n'% load[i])
        file2.write(l)
    file2.flush()
    file2.close
    return

def write_qcal(obsid):
    qcal=['#!/bin/bash -l']
    qcal.append('#SBATCH --nodes=25')
    qcal.append('#SBATCH --ntasks-per-node=1')
    qcal.append('#SBATCH --time=1:00:00')
    qcal.append('#SBATCH --partition=gpuq')
    qcal.append('#SBATCH --account=mwasci')
    qcal.append('#SBATCH --export=NONE')
    qcal.append('#SBATCH --output="%s.log"'% obsid)
    qcal.append('#SBATCH --open-mode=append')
    qcal.append('m -f cal.log')
    qcal.append('./%scal.py'% obsid)
    #qcal.append('/home/ibrown/bin/bpsummary.py --obsid=%s'% obsid)
    
    #write qcal.sh
    file1 = open('%sqcal.sh'% obsid,'w')
    for i in range (0, len(qcal)):
        q='%s\n'% qcal[i]
        file1.write(q) 
    file1.flush()
    file1.close
    return

def write_cal(field_size,fscrunch,img_cad, obsid):
    cal=['#!/usr/bin/env python']
    cal.append('from mwa_pipe import *\n')
    cal.append('import timeit')
    cal.append('import os')
    cal.append('start = timeit.default_timer()')
    cal.append('sys.stdout.write("***********************************\\n")')
    cal.append('sys.stdout.write("Begin calibrating observation data \\n")')
    cal.append('sys.stdout.write("***********************************\\n")')
    cal.append('sys.stdout.flush()')
    cal.append('rts = MWAPipe("cal")\n')
    cal.append('rts.field_size = %d'% field_size)
    cal.append('rts.update_cal = True')
    cal.append('rts.ncalibrators = 1')
    cal.append('rts.template_base = "rts_"')
    cal.append('rts.fscrunch = %d'% fscrunch)
    cal.append('rts.cat_extent = 20.0')
    cal.append('rts.cal_cadence = 64')
    cal.append('rts.img_cadence = %d'% img_cad)
    cal.append('rts.npeel = 0')
    cal.append('rts.niono_cal = 0')
    cal.append('rts.cal_srcs = 30')
    cal.append('rts.rts_bin = "rts_gpu"')
    #cal.append('obsid = get_local_obs_list(%s)'% obsid)
    #cal.append('for obsid in obsid:')
    cal.append('# Clean out any old files')
    cal.append('rts.clean("%s")'% obsid)
    cal.append('# Calibrate the obsid')
    cal.append('rts.calibrate("%s")'% obsid)
    cal.append('os.system("/home/ibrown/bin/bpsummary.py --obsid=%s")'% obsid)
    cal.append('sys.stdout.write("***********************************\\n")')
    cal.append('sys.stdout.write("Finished calibrating observation data \\n")')
    cal.append('sys.stdout.write("***********************************\\n")')
    cal.append('print("total time to calibrate data:%f")%(timeit.default_timer()-start)')
    
    #write cal.py
    file2 = open('%scal.py'% obsid,'w')
    for i in range (0,len(cal)):
        c=('%s\n'% cal[i])
        file2.write(c)
    file2.flush()
    file2.close
    return

def write_qimage (weighting,obsid, email):
    qimage=['#!/bin/bash -l']
    qimage.append('#SBATCH --nodes=25')
    qimage.append('#SBATCH --ntasks-per-node=1')
    qimage.append('#SBATCH --time=1:00:00')
    
    if weighting=='natural':
        qimage.append('#SBATCH --partition=gpuq')
    else:
        qimage.append('#SBATCH --partition=workq')
        
    qimage.append('#SBATCH --account=mwasci')
    qimage.append('#SBATCH --export=NONE')
    if email != 'none':
        qimage.append('#SBATCH --mail-type=END')
        qimage.append('#SBATCH --mail-user=%s'% email)
    qimage.append('#SBATCH --output="%s.log"'% obsid)
    qimage.append('#SBATCH --open-mode=append')
    qimage.append('rm -f image.log')
    qimage.append('./%simage.py'% obsid)
    
    
    #write qimage.sh
    file3 = open('%sqimage.sh'% obsid,'w')
    for i in range (0,len(qimage)):
        qimg=('%s\n'%qimage[i])
        file3.write(qimg)
    file3.flush()
    file3.close
    return

def write_image(obsid, weighting,robustness,field_size,fscrunch,img_cad,path):
    image=['#!/usr/bin/env python']
    image.append('from mwa_pipe import *\n')
    image.append('import timeit')
    image.append('rts = MWAPipe("image")\n')
    image.append('#')
    image.append('rts.weighting = "%s"'% weighting)
    
    if weighting=='robust':
        image.append('rts.robustness = %s'% robustness)
    
    image.append('rts.field_size = %s'% field_size)
    image.append('rts.ncalibrators = 1')
    image.append('rts.template_base = "rts_"')
    image.append('rts.fscrunch = %s'% fscrunch)
    image.append('rts.cat_extent = 20.0')
    image.append('rts.cal_cadence = 64')
    image.append('rts.img_cadence = %s'% img_cad)
    image.append('rts.cal_srcs = 30')
    image.append('rts.do_accumulate = False')
    image.append('#')
    image.append('# Test set for transient')
    image.append('rts.min_baseline = 50.0')
    image.append('rts.make_psf = False')
    image.append('rts.cat_extent = 6.0')
    image.append('rts.npeel = 0')
    image.append('rts.niono_cal = 0')
    image.append('rts.img_srcs = 100')
    image.append('rts.update_cal = False')
    image.append('#')
    
    if weighting=='natural':
        image.append('rts.rts_bin = "rts_gpu"')
    else:
        image.append('rts.rts_bin = "rts_cpu"')
    image.append('rts.img_cat = None')
    image.append('# Image the obsid using in-field calibration performed earlier')
    image.append('start = timeit.default_timer()')
    image.append('sys.stdout.write("***********************************\\n")')
    image.append('sys.stdout.write("Begin imaging observation data \\n")')
    image.append('sys.stdout.write("***********************************\\n")')
    image.append('sys.stdout.flush()')
    image.append('rts.image(%s, cal_id = None, dest_image_path = "%s")'% (obsid, path))
    image.append('sys.stdout.write("***********************************\\n")')
    image.append('sys.stdout.write("Finished imaging observation data \\n")')
    image.append('sys.stdout.write("***********************************\\n")')
    image.append('print("total time to image data:%f")%(timeit.default_timer()-start)')
    
    #write image.py
    file4 = open('%simage.py'% obsid,'w')
    for i in range (0,len(image)):
        img=('%s\n'% image[i])
        file4.write(img)
    file4.flush()
    file4.close
    return

def write_qchkcal (obsid,plot):
    qchkcal=['#!/bin/bash -l']
    qchkcal.append('#SBATCH --nodes=1')
    qchkcal.append('#SBATCH --ntasks=1')
    qchkcal.append('#SBATCH --ntasks-per-node=1')
    if plot=='true':
        qchkcal.append('#SBATCH --time=1:00:00')
    else:
        qchkcal.append('#SBATCH --time=0:15:00')
    qchkcal.append('#SBATCH --partition=gpuq')
    qchkcal.append('#SBATCH --account=mwasci')
    qchkcal.append('#SBATCH --export=NONE')
    qchkcal.append('#')
    qchkcal.append('#')
    qchkcal.append('/home/ibrown/bin/bpsummary.py --obslist=%s --plot=%s'% (obsid,plot))
    
    #write qchkcal
    file = open('qchkcal.sh','w')
    for i in range (0,len(qchkcal)):
        chk=('%s\n'% qchkcal[i])
        file.write(chk)
    file.flush()
    file.close
    return

def write_rts_cal():
    rtscal=['//---------//\n']
    rtscal.append('FscrunchChan=4\n')
    rtscal.append('//---------//\n')
    rtscal.append('SubBandIDs=\n')
    rtscal.append('StorePixelMatrices=0\n')
    rtscal.append('// do extra iterations with in each time step')
    rtscal.append('//IterateOverAllGains=0\n')
    rtscal.append('MaxFrequency=180')
    rtscal.append('ImageOversampling=3\n')
    rtscal.append('applyDIcalibration=1')
    rtscal.append('doMWArxCorrections=1')
    rtscal.append('doRFIflagging=1')
    rtscal.append('useFastPrimaryBeamModels=1\n')
    rtscal.append('CorrDumpTime=0.0')
    rtscal.append('CorrDumpsPerCadence=0v')
    rtscal.append('NumberOfIntegrationBins=0')
    rtscal.append('NumberOfIterations=\n')
    rtscal.append('StartProcessingAt=0')
    rtscal.append('StartIntegrationAt=0\n')
    rtscal.append('// -------------------------------------------------------------------------- //')
    rtscal.append('// Need to specify extra info when not using uvfits files\n')
    rtscal.append('// in correlator mode, Base File name is used to find correlator files.') 
    rtscal.append('BaseFilename=*_gpubox\n')
    rtscal.append('doRawDataCorrections=1\n')
    rtscal.append('ReadGpuboxDirect=1')
    rtscal.append('//UseCorrelatorInput=1')
    rtscal.append('UsePacketInput=0')
    rtscal.append('UseThreadedVI=0\n')
    rtscal.append('ArrayFile=array_file.txt')
    rtscal.append('ArrayNumberOfStations=128\n')
    rtscal.append('ChannelBandwidth=0.04')
    rtscal.append('NumberOfChannels=32\n')
    rtscal.append('ArrayPositionLat=-26.70331940')
    rtscal.append('ArrayPositionLong=116.67081524\n')
    rtscal.append('//time is needed to set lst, even if ObservationTimeBase is set.')
    rtscal.append('ObservationTimeBase=2456519.20083\n')
    rtscal.append('// -------------------------------------------------------------------------- //\n')
    rtscal.append('// --- observing stuff --- //\n')
    rtscal.append('ReadAllFromSingleFile=\n')
    rtscal.append('ObservationFrequencyBase=167.055\n')
    rtscal.append('ObservationPointCentreHA=1.031')
    rtscal.append('ObservationPointCentreDec=-25.93\n')
    rtscal.append('ObservationImageCentreRA=')
    rtscal.append('ObservationImageCentreDec=\n')
    rtscal.append('// --- weighting stuff --- //\n')
    rtscal.append('calBaselineMin=20.0')
    rtscal.append('calShortBaselineTaper=50.0\n')
    rtscal.append('// --- calibration stuff --- //\n')
    rtscal.append('DoCalibration=\n')
    rtscal.append('SourceCatalogueFile=catalogue.txt')
    rtscal.append('NumberOfCalibrators=1\n')
    
    #write rts_cal.in
    file = open('rts_cal.in','w')
    for i in range (0,len(rtscal)):
        rcal=('%s\n'% rtscal[i])
        file.write(rcal)
    file.flush()
    file.close
    return

def write_rts_img(maxfreq=180,oversample=5):
    rtsimg=['//---------//\n']
    rtsimg.append('//ImagePSF=1')
    rtsimg.append('FscrunchChan=4\n')
    rtsimg.append('//---------//\n')
    rtsimg.append('SubBandIDs=\n')
    rtsimg.append('StorePixelMatrices=1\n')
    rtsimg.append('// do extra iterations with in each time step')
    rtsimg.append('//IterateOverAllGains=0\n')
    rtsimg.append('MaxFrequency=%d'%maxfreq)
    rtsimg.append('ImageOversampling=%d\n'% oversample)
    rtsimg.append('applyDIcalibration=1')
    rtsimg.append('doMWArxCorrections=1')
    rtsimg.append('doRFIflagging=0')
    rtsimg.append('useFastPrimaryBeamModels=1\n')
    rtsimg.append('CorrDumpTime=0.0\n')
    rtsimg.append('CorrDumpsPerCadence=0')
    rtsimg.append('NumberOfIntegrationBins=0')
    rtsimg.append('NumberOfIterations=\n')
    rtsimg.append('StartProcessingAt=0')
    rtsimg.append('StartIntegrationAt=0\n')
    rtsimg.append('// -------------------------------------------------------------------------- //')
    rtsimg.append('// Need to specify extra info when not using uvfits files\n')
    rtsimg.append('// in correlator mode, Base File name is used to find correlator files.')
    rtsimg.append('ReadAllFromSingleFile=')
    rtsimg.append('BaseFilename=*_gpubox\n')
    rtsimg.append('doRawDataCorrections=1\n')
    rtsimg.append('ReadGpuboxDirect=1')
    rtsimg.append('//UseCorrelatorInput=1')
    rtsimg.append('UsePacketInput=0')
    rtsimg.append('UseThreadedVI=1\n')
    rtsimg.append('ArrayFile=array_file.txt')
    rtsimg.append('ArrayNumberOfStations=128\n')
    rtsimg.append('ChannelBandwidth=0.04')
    rtsimg.append('NumberOfChannels=32\n')
    rtsimg.append('ArrayPositionLat=-26.70331940')
    rtsimg.append('ArrayPositionLong=116.67081524\n')
    rtsimg.append('//time is needed to set lst, even if ObservationTimeBase is set.')
    rtsimg.append('ObservationTimeBase=2456519.20083\n')
    rtsimg.append('// -------------------------------------------------------------------------- //\n')
    rtsimg.append('// --- observing stuff --- //\n')
    rtsimg.append('ObservationFrequencyBase=167.055\n')
    rtsimg.append('ObservationPointCentreHA=1.031')
    rtsimg.append('ObservationPointCentreDec=-25.93\n')
    rtsimg.append('ObservationImageCentreRA=')
    rtsimg.append('ObservationImageCentreDec=\n')
    rtsimg.append('// --- weighting stuff --- //\n')
    rtsimg.append('calBaselineMin=20.0')
    rtsimg.append('calShortBaselineTaper=50.0\n')
    rtsimg.append('// --- calibration stuff --- //\n')
    rtsimg.append('DoCalibration=\n')
    rtsimg.append('SourceCatalogueFile=catalogue_imaging.txt\n')
    rtsimg.append('NumberOfCalibrators=1')
    rtsimg.append('NumberOfSourcesToPeel=0')
    rtsimg.append('NumberOfIonoCalibrators=0')
    rtsimg.append('UpdateCalibratorAmplitudes=0\n')
    
    #write rts_img.in
    file = open('rts_img.in','w')
    for i in range (0,len(rtsimg)):
        rimg=('%s\n'% rtsimg[i])
        file.write(rimg)
    file.flush()
    file.close
    return

def write_rts_uv():
    rtsuv=['//---------//\n']
    rtsuv.append('//ImagePSF=1')
    rtsuv.append('FscrunchChan=4\n')
    rtsuv.append('//---------//\n')
    rtsuv.append('SubBandIDs=\n')
    rtsuv.append('//StorePixelMatrices=1\n')
    rtsuv.append('// do extra iterations with in each time step')
    rtsuv.append('//IterateOverAllGains=0\n')
    rtsuv.append('MaxFrequency=200')
    rtsuv.append('ImageOversampling=5')
    rtsuv.append('writeVisToUVFITS=1\n')
    rtsuv.append('applyDIcalibration=1')
    rtsuv.append('doMWArxCorrections=1')
    rtsuv.append('doRFIflagging=0')
    rtsuv.append('useFastPrimaryBeamModels=1\n')
    rtsuv.append('CorrDumpTime=0.0\n')
    rtsuv.append('CorrDumpsPerCadence=0')
    rtsuv.append('NumberOfIntegrationBins=0')
    rtsuv.append('NumberOfIterations=\n')
    rtsuv.append('StartProcessingAt=0')
    rtsuv.append('StartIntegrationAt=0\n')
    rtsuv.append('// -------------------------------------------------------------------------- //')
    rtsuv.append('// Need to specify extra info when not using uvfits files\n')
    rtsuv.append('// in correlator mode, Base File name is used to find correlator files. ')
    rtsuv.append('ReadAllFromSingleFile=')
    rtsuv.append('BaseFilename=*_gpubox\n')
    rtsuv.append('doRawDataCorrections=1\n')
    rtsuv.append('ReadGpuboxDirect=1')
    rtsuv.append('//UseCorrelatorInput=1')
    rtsuv.append('UsePacketInput=0')
    rtsuv.append('UseThreadedVI=1\n')
    rtsuv.append('ArrayFile=array_file.txt')
    rtsuv.append('ArrayNumberOfStations=128\n')
    rtsuv.append('ChannelBandwidth=0.04')
    rtsuv.append('NumberOfChannels=32\n')
    rtsuv.append('ArrayPositionLat=-26.70331940')
    rtsuv.append('ArrayPositionLong=116.67081524\n')
    rtsuv.append('//time is needed to set lst, even if ObservationTimeBase is set.')
    rtsuv.append('ObservationTimeBase=2456519.20083\n')
    rtsuv.append('// -------------------------------------------------------------------------- //\n')
    rtsuv.append('// --- observing stuff --- //\n')
    rtsuv.append('ObservationFrequencyBase=167.055\n')
    rtsuv.append('ObservationPointCentreHA=1.031')
    rtsuv.append('ObservationPointCentreDec=-25.93\n')
    rtsuv.append('ObservationImageCentreRA=')
    rtsuv.append('ObservationImageCentreDec=\n')
    rtsuv.append('// --- weighting stuff --- //\n')
    rtsuv.append('calBaselineMin=20.0')
    rtsuv.append('calShortBaselineTaper=50.0\n')
    rtsuv.append('//imgLongBaselineTaper = 3000.0')
    rtsuv.append('imgBaselineMin = 100.0\n')
    rtsuv.append('// --- calibration stuff --- //\n')
    rtsuv.append('DoCalibration=\n')
    rtsuv.append('SourceCatalogueFile=catalogue_imaging.txt\n')
    rtsuv.append('NumberOfCalibrators=1')
    rtsuv.append('NumberOfSourcesToPeel=0')
    rtsuv.append('NumberOfIonoCalibrators=0')
    rtsuv.append('UpdateCalibratorAmplitudes=0\n')
    rtsuv.append('// --- imaging stuff --- //\n')
    rtsuv.append('// Turn on imaging\n')
    rtsuv.append('//MakeImage=\n')
    rtsuv.append('//FieldOfViewDegrees=5\n')
    rtsuv.append('//GridDataWithMethod=2\n')
    rtsuv.append('//MakeStokesSnapshots=')
    rtsuv.append('//MakeWeightedSnapshots=\n')
    rtsuv.append('//RegridMethod=2')
    rtsuv.append('//DoRegriddingWithProjection=2048\n')
    
        #write rts_img.in
    file = open('rts_uv.in','w')
    for i in range (0,len(rtsuv)):
        ruv=('%s\n'% rtsuv[i])
        file.write(ruv)
    file.flush()
    file.close
    return
