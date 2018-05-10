#!/usr/bin/python
#TODO:
#pick dynamic max value for re-normalization

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
from biomarker_labeling import labelBiomarkers
from processSegmentation import 

#Take plate folder as input 



#inputPath = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/insets_TxB';    
#outputPath = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/output';
#numCycles = 11; %Is this basically number of cycles? 2:cycle num
#row = ['B':'D'];   
#col = [2:5];
#SaveFig = 1;

#Inputting folder for now: assume mix of rows and cols.
#will have to think more about how to handle subfolders, plates, any other organization. 
#Tenuous but potentially useful assumption: first character is row, 2/3 is col 
#can potentially search based on rows/cols in input file. 



def resegment(imageFile):

    print('/n/scratch2/rps21/ResegmentedPlates/%s' % platename, nuclearstack, rows, cols)
    subprocess.call(["/home/rps21/pCycIF_Segmentation/run_runSegmentationO2.sh /home/rps21/MCR/v901/ externalParams.txt"])

#    print('done segmenting %s' % plate)



#with open('../externalParams.txt') as paramFile:
#    paramLines = paramFile.readlines()
#for line in paramLines:
#    if 'inputPath' in line and '=' in line:
#        inputLine = line
#folderToSegment = inputLine.split('=')[1].split('%')[0].strip()
#folderToSegment = folderToSegment.strip('\'')
#for fn in sorted(glob.glob(folderToSegment+'/*.tif')):
#    resegment(fn)

#fill biomarkers
#run normalization
#combine fields.
#organize files. 

resegment_plate(platename)


