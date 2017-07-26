#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 11:16:50 2017

@author: ianbrown
"""
'''
fscrunch=32
field_size=30.0
obsid='1165925976'
img_cad= 96
path='test%s'% obsid
weighting='uniform'
robustness=(-1)
plot='true'
'''
import commands, os, time, sys
from optparse import OptionParser
import rtswrite as rtsw

 
parser = OptionParser()#create parser for program options

parser.add_option("-o", "--obsid",type='string',action="store", dest="obsid",help='obsid to be used')
#parser.add_option("-p","--plot", type='string',action='store', dest="plot", help='plot calibration data', default='false')
parser.add_option("-i","--imgcad",type='float', action='store',dest='img_cad',help='image cadence', default=96)
parser.add_option("-w","--weighting",type='string',action='store',dest='weighting',help='type of weighting to use',default='uniform')
parser.add_option('-r','--robustness',type='float',action='store',dest='robustness',help='value for robust weighting',default=-1.0)
parser.add_option('-p','--path',type='string', action= 'store', dest='path',help='path to creat fits files in', default='default')
parser.add_option('-s','--size',type='float',action='store',dest='field_size',help='field of view in degrees',default=5.0)
parser.add_option('--fscrunch',type='int', action='store', dest='fscrunch', help='number of channels to average together',default=32)
parser.add_option('-e','--email',type='string',action='store',dest='email', help='email adress', default='none')
parser.add_option('--oversample',type='int',action='store',dest='oversample', help='number of pixels per beam', default=5)
parser.add_option('-m','--maxfreq',type='int',action='store',dest='maxfreq', help='maximum frequency to set beam size', default=180)

options, args = parser.parse_args()
#main_path=os.getcwd() # set main directory path

#write files for rts execution
rtsw.write_qload(options.obsid,options.email)
rtsw.write_qcal(options.obsid)
rtsw.write_qimage(options.weighting,options.obsid, options.email)
#write_qchkcal(options.obsid,options.plot)
rtsw.write_load(options.obsid)
rtsw.write_cal(options.field_size,options.fscrunch,options.img_cad, options.obsid)
rtsw.write_image(options.obsid, options.weighting,options.robustness,options.field_size,options.fscrunch,options.img_cad,options.path)
rtsw.write_rts_img(options.maxfreq,options.oversample)
rtsw.write_rts_cal()
rtsw.write_rts_uv()

os.system('chmod u+x ./*.py')#make .py programs executable

# submit the load job
cmd = "sbatch %sqload.sh"%(options.obsid)
print ("Submitting load with command: %s" % cmd)
status, jobnum = commands.getstatusoutput(cmd)
if (status == 0 ):
    print ("%s" % jobnum)
    jnum=jobnum.split(' ')
    ql_id=(jnum[3])[:7]
    
else:
    print ("Error submitting load")
    sys.exit('failed')

#wait for files to load
time.sleep(120)
#check for completion of load process
status='RUNNING'
while status != 'COMPLETED':
    time.sleep(60)
    cmd = "sacct --cluster=zeus -j %s"% ql_id
    print ("Submitting sacct with command: %s" % cmd)
    exitcode, jobnfo = commands.getstatusoutput(cmd)
    jnfo=jobnfo.split() 
    status=jnfo[19]
    if status == 'FAILED':
        sys.exit(status)
            
#submit cal job
cmd = "sbatch %sqcal.sh"%(options.obsid)
print ("Submitting cal with command: %s" % cmd)
status, jobnum = commands.getstatusoutput(cmd)
if (status == 0 ):
    print ("%s" % jobnum)
    jnum=jobnum.split(' ')
    cal_id=(jnum[3])[:7]
    #submit check calibration    
    #cmd = "sbatch --dependency=afterok:%s qchkcal.sh"% cal_id
    #print ("Submitting chk_cal with command: %s" % cmd)
    #status, jobnum = commands.getstatusoutput(cmd)
    #if (status == 0 ):
     #   print ("%s" % jobnum)
     #   jnum=jobnum.split(' ')
     #   chkcal_id=(jnum[3])[:7]
     #submit image job
    cmd = "sbatch --dependency=afterok:%s %sqimage.sh"% (cal_id,(options.obsid))
    print ("Submitting image with command: %s" % cmd)
    status, jobnum = commands.getstatusoutput(cmd)
    if (status == 0 ):
        print ("%s" % jobnum)
        jnum=jobnum.split(' ')
        img_id=(jnum[3])[:7]
        
    else:
        print ("Error submitting qimage")    
    #else:
     #   print ("Error submitting qchkcal")
else:
    print ("Error submitting qcal")