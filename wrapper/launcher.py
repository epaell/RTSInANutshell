#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:51:33 2017
RTS wrapper launch
@author: ianbrown
"""

import commands,sys
from rtswrite import write_launch
fscrunch=1
field_size=8.0
obslist=['1165925976']
img_cad= 96
path='uvceti'
weighting='robust'
robustness=(-1)
email='ianbrown@uwm.edu'
oversample=5
maxfreq=180


#main launch job string for each job id
for obsid in obslist:
    write_launch(obsid,fscrunch,img_cad,path,weighting,robustness,field_size,email,oversample,maxfreq)
    # submit the load job
    cmd = "sbatch %srts.sh"%(obsid)
    print ("Submitting load with command: %s" % cmd)
    status, jobnum = commands.getstatusoutput(cmd)
    if (status == 0 ):
        print ("%s" % jobnum)
                
    else:
        print ("Error submitting %sr"% obsid)
        sys.exit('failed')
