#!/usr/bin/python

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
import yaml
import os
from processSegmentation import labelBiomarkers
from processSegmentation import normalizePlate



subprocess.call(['./run_runSegmentation_Caitlin.sh','/opt/mcr/v94/','config.yaml'])

with open('config.yaml') as stream:
    options = yaml.load(stream)
outputPath = options['outputPath']
numCycles = options['numCycles']
markerList = options['markerList']
bleached = options['bleachImaged']

newFolder = outputPath+'/processedOutput/'
if not os.path.exists(newFolder):
    os.makedirs(newFolder)

for fn in sorted(glob.glob(outputPath+'/*.txt')):
    #fill biomarkers
    #df = labelBiomarkers(fn,numCycles,markerList,bleached)
    filename_out = newFolder + 'name' + fn.split('/')[-1]
    df.to_csv(filename_out)
    #clean and log-transform data
    df = normalizePlate(df)   #create new output subdirectory?
    #Save file
    filename_out = newFolder + 'norm' + fn.split('/')[-1]
    df.to_csv(filename_out)


#add nuclei thresholding?

