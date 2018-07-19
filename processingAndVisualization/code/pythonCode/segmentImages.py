#!/usr/bin/python

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
import yaml
import os
from biomarker_labeling import labelBiomarkers
from processSegmentation import normalizePlate



subprocess.call(['./run_runSegmentation.sh','/opt/mcr/v94/','config.yaml'])

with open('config.yaml') as stream:
    options = yaml.load(stream)
outputPath = options['outputPath']
numCycles = options['numCycles']


newFolder = outputPath+'/processedOutput/'
if not os.path.exists(newFolder):
    os.makedirs(newFolder)

for fn in sorted(glob.glob(outputPath+'/*.txt')):
    #fill biomarkers
    #df = labelBiomarkers(fn,numCycles,markerList=None,bleached=1)
    #now in normalize plate. 
    #run normalization
    df=normalizePlate(fn,numCycles,markerList=None)   #create new output subdirectory?
    #Save file
    filename_out = newFolder + fn.split('/')[-1]
    df.to_csv(filename_out)




#add nuclei thresholding?
#markerList = ['DNA-1','DNA-2','DNA-3','DNA-4','DNA-5','MLC2v','pS6-240','p-Src','c-Jun','CleavedCaspase3','pan14-3-3','pRb','Oct4a','B-actin','B-tubulin','pAKT-S473','pS6-235','p-mTor','NFkB','mTor'] #
