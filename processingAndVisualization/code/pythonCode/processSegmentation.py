#!/usr/bin/python

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
import os
from copy import deepcopy


#Need checks for length based on cycle number, check new col length vs old. 
def labelBiomarkers(fn,numCycles,markerList=None,bleached=1):
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
            columnLabels.append(markerList[i]+'_Nuc')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_Nuc')
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_Cyto')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_Cyto')


    #If list is provided, use for labeling
    else:
        channelList = ['DAPI','FITC','Cy3','Cy5']
        for i in range(1,numCycles+1):
#            for j in range(0,numCycles):
#                dapiList.append('DNA-%s' % str(j+1) + '_Nuc')
            for channel in channelList:
                label ='%s-%s' % (channel,i)
                columnLabels.append(label + '_Nuc')
        for i in range(1,numCycles+1):
#            for j in range(0,numCycles):
#                dapiList.append('DNA-%s' % str(j+1) + '_Cyto')
            for channel in channelList:
                label ='%s-%s' % (channel,i)
                columnLabels.append(label + '_Cyto')
    #missing and bleach 

    #For now hard code, may want to read these from original file, or get input from user
    extraColumns = ['NucleusArea','CytoplasmArea','CellPosition_X','CellPosition_Y']
    columnLabels = columnLabels + extraColumns

    originalDF = pd.read_table(fn)
    outputDF = deepcopy(originalDF)
    outputDF.columns = columnLabels #Here. column labels must be wrong. 
#    newFN = fn.split('.')[0]+'.csv'
#    outputDF.to_csv(newFN)

    return outputDF



def normalizePlate(df):
    #Eventually include histogram stretching
    df = df.replace('', np.nan)
    df = df.replace(0, np.nan)  #Could this have any negative consequences?
    df = df.replace(float('-inf'), np.nan)
    df = df.replace(np.inf, np.nan)
#    df = df.dropna() #default options should be appropriate

    #Include optional list of non intesnity columns
    #For now, know which columns exist from typical segmentation output 
    nonIntCols = ['NucleusArea','CytoplasmArea','CellPosition_X','CellPosition_Y']
#    optDF = pd.DataFrame() 
    for col in df.columns:
        if col not in nonIntCols:
            df[col] = df[col].apply(np.log2)



#    pd.concat([df,optDF], axis=1)
    #optionally save at each step?:
#        filename_out = normalizedOutputFolder + 'normalized' + filename.split('/')[-1]
#        df.to_csv(filename_out)

    return df


#TODO:
#Histograms stretching
                #Eventually add back in automated histogram stretching
                #Need automated way to find 'newmax' for each plate 
#                newmin = min(df[col])
#                newmax = 15000  #Really need a better way of setting this
#                oldmax = max(df[col])
#                if (oldmax - newmin) == 0:
#                    print(fn)
#                    df[col] = np.NaN
#                else:
#                    df[col] = (df[col] - newmin)*((newmax-newmin)/(oldmax - newmin))+newmin
 

#Option dict input for drug, dose, time
            #add dose and drug columns
#            doseDrug = filename.split('_')[1]
#            row = doseDrug[0] #Need to convert string to num into final csv
#            if plate == 'ProcessedDataP6':
#                dose_dict = {('B','proteome1'):10,('D','proteome2'):1,('F','proteome3'):3,('C','proteome1'):3,('E','proteome2'):10,('G','proteome3'):1}
#            else:
#                dose_dict = {('B','signaling1'):10,('E','signaling2'):10,('C','signaling1'):3,('F','signaling2'):3,('B','proteome1'):3,('D','proteome2'):3,('F','proteome3'):3,('D','signaling1'):1,('G','signaling2'):1,('C','proteome1'):1,('E','proteome2'):1,('G','proteome3'):1}
#            dose = dose_dict[(row,ab_set)]
#            df['Dose'] = pd.Series([dose]*len(df.index),index = df.index)
#    
#            drug = doseDrug[1:]
#            df['Drug'] = pd.Series([drug]*len(df.index),index = df.index)


#multiple output folders for various analyses 

#def setOutput(paramsFN):
#    with open('externalParams.txt') as paramFile:   #File location is an issue
#        paramLines = paramFile.readlines()
#    for line in paramLines:
#        if 'outputPath' in line and '=' in line:
#            outputLine = line
#    outputFolder = outputLine.split('=')[1].split('%')[0].strip()
#    outputFolder = outputFolder.strip('\'')

#    normalizedOutput = outputFolder+'/normalized/'
#    if not os.path.exists(normalizedOutput):
#        os.makedirs(normalizedOutput)
#    print(normalizedOutput)

#    folders = [outputFolder,normalizedOutput]    

#    return folders




