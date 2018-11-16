#!/usr/bin/python

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
import os
from copy import deepcopy


#TODO: Need checks for new col length vs old, handle error


def labelBiomarkers(fn,numCycles,markerList=None,bleached=0):
    columnLabels = []
    numMarkers = numCycles*4
    dapiList = []
    bleached=0
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
#    try:
#        outputDF.columns = columnLabels 
#    except ValueError:
    if len(outputDF.columns) != len(columnLabels):
        raise ValueError('Incorrect number of markers or cycles given. Expected %s markers, received %s' % (len(outputDF.columns),len(columnLabels)))
    else:
        outputDF.columns = columnLabels 
    return outputDF





#TODO:
#Histograms stretching
                #Need automated way to find 'newmax' for each plate 
#                newmin = min(df[col])
#                newmax = 15000  #Need to automate finding this value
#                oldmax = max(df[col])
#                if (oldmax - newmin) == 0:
#                    print(fn)
#                    df[col] = np.NaN
#                else:
#                    df[col] = (df[col] - newmin)*((newmax-newmin)/(oldmax - newmin))+newmin
 






