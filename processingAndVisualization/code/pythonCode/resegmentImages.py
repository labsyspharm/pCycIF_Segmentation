#!/usr/bin/python

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
from biomarker_labeling import labelBiomarkers
from processSegmentation import normalizePlate
from processSegmentation import setOutput

#    subprocess.call(["/home/rps21/Cardiomyocite/segmentation/bobbySegmentCardioMyo/run_bobbySegmentCardiomyo_orchestra_addbackdapi.sh", "/home/rps21/MCR/v901/", '/n/scratch2/rps21/ResegmentedPlates/%s' % platename, nuclearstack, rows,

subprocess.call(['/home/rps21/pCycIF_Segmentation/run_runSegmentationO2.sh','/home/rps21/MCR/v901/','externalParams.txt'])
outputFolders = setOutput('externalParams.txt')
#Need to read these inputs more effectively. 
numCycles = 6 
markerList = ['Rb-MAVS','LAMP2','M-IRF3','p-TBK1','NFAT-C1','NFkB','TMEM137','IRF5','COX4','STAT6','IRF7','IRF1','STAT5b','LC3A/B','STAT5a','Stat3','PKM2','p-mTOR']
normalizePlate(outputFolders,markerList,numCycles)





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



#def resegment(imageFile):

#    print('/n/scratch2/rps21/ResegmentedPlates/%s' % platename, nuclearstack, rows, cols)
#    subprocess.call(["/home/rps21/pCycIF_Segmentation/run_runSegmentationO2.sh /home/rps21/MCR/v901/ externalParams.txt"])




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




