#!/usr/bin/python

import pandas as pd
import glob
import sys 
import subprocess
import yaml
import os
import numpy as np
import os
from copy import deepcopy


def labelBiomarkers(fn,numCycles,markerList=None,bleached=0):
    columnLabels = []
    numMarkers = numCycles*4
    dapiList = []
    
    #If no list of markers is provided, create generic naming
    if markerList:
        #add dapi labels to marker list. May need to handle this being optional 
        for i in range(0,numCycles):
            dapiList.append('DNA-%s' % str(i+1))
        markerList = dapiList+markerList
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_MeanIntensity_Nuc')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_MeanIntensity_Nuc')
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_MeanIntensity_Cyto')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_MeanIntensity_Cyto')
        #Now including median as well, repeat loops for now. May want to do something based off initial column names
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_MedianIntensity_Nuc')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_MedianIntensity_Nuc')
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_MedianIntensity_Nuc')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_MedianIntensity_Nuc')




    #If list is provided, use for labeling
    else:
        channelList = ['DAPI','FITC','Cy3','Cy5']
        for i in range(1,numCycles+1):
            for channel in channelList:
                label ='%s-%s' % (channel,i)
                columnLabels.append(label + '_MeanIntensity_Nuc')
                if bleached == 1 and (i+1) % numCycles != 0:   
                    columnLabels.append('Bleached-%s' % markerList[i] + '_MeanIntensity_Nuc')
        for i in range(1,numCycles+1):
            for channel in channelList:
                label ='%s-%s' % (channel,i)
                columnLabels.append(label + '_MeanIntensity_Cyto')
                if bleached == 1 and (i+1) % numCycles != 0:   
                    columnLabels.append('Bleached-%s' % markerList[i] + '_MeanIntensity_Cyto')
        #Now including median as well, repeat loops for now. May want to do something based off initial column names
        for i in range(1,numCycles+1):
            for channel in channelList:
                label ='%s-%s' % (channel,i)
                columnLabels.append(label + '_MedianIntensity_Nuc')
                if bleached == 1 and (i+1) % numCycles != 0:   
                    columnLabels.append('Bleached-%s' % markerList[i] + '_MedianIntensity_Nuc')
        for i in range(1,numCycles+1):
            for channel in channelList:
                label ='%s-%s' % (channel,i)
                columnLabels.append(label + '_MedianIntensity_Cyto')
                if bleached == 1 and (i+1) % numCycles != 0:   
                    columnLabels.append('Bleached-%s' % markerList[i] + '_MedianIntensity_Cyto')



    #For now hard code, may want to read these from original file, or get input from user
    extraColumns = ['NucleusArea','CytoplasmArea','CellPosition_X','CellPosition_Y']
    columnLabels = columnLabels + extraColumns


    originalDF = pd.read_table(fn)
    outputDF = deepcopy(originalDF)

    if len(outputDF.columns) != len(columnLabels):
        raise ValueError('Incorrect number of markers or cycles given. Expected %s markers, received %s' % (len(outputDF.columns),len(columnLabels)))
    else:
        outputDF.columns = columnLabels 
    return outputDF

#segment cells
subprocess.call(['/segmentation/run_runSegmentation.sh','/opt/mcr/v94/','/config/cycif_segmentation.yml'])

with open('/config/cycif_segmentation.yml') as stream:
    options = yaml.load(stream)
outputPath = '/output'
numCycles = options['numCycles']
markerList = options['markerList']
bleached = options['bleachImaged']

for fn in sorted(glob.glob(outputPath+'/*.txt')):
    df = labelBiomarkers(fn,numCycles,markerList,bleached)
    df.to_csv(fn)
