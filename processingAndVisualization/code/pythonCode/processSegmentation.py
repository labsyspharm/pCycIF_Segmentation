#!/usr/bin/python

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
import os
from biomarker_labeling import labelBiomarkers

#This gets file names
#To start want to we write everything with same first three characters (row/col) to one df/file 
#get from fn.split('/')[-1].split('_')[0]

def setOutput(paramsFN):
    with open('externalParams.txt') as paramFile:   #File location is an issue
        paramLines = paramFile.readlines()
    for line in paramLines:
        if 'outputPath' in line and '=' in line:
            outputLine = line
    outputFolder = outputLine.split('=')[1].split('%')[0].strip()
    outputFolder = outputFolder.strip('\'')

    normalizedOutput = outputFolder+'/normalized/'
    if not os.path.exists(normalizedOutput):
        os.makedirs(normalizedOutput)
    print(normalizedOutput)

    folders = [outputFolder,normalizedOutput]    

    return folders



def normalizePlate(outputFolders,markerList,numCycles):

    outputFolder = outputFolders[0]
    normalizedOutputFolder = outputFolders[1]

    #mkdir for normalized output      
    for fn in sorted(glob.glob(outputFolder+'/_*.txt')):
#        df = pd.read_table(fn)
        df = labelBiomarkers(fn,numCycles,markerList,1)
        #Normalize data
        for col in df.columns:

            #Include optional list of non intesnity columns
            #Will save them in separate df, exclude from normalization, add back in later.
            #May be good to automate searching for this type of column, but for now can feed a list.
            nonIntCols = ['NucleusArea','CytoplasmArea','CellPosition_X','CellPosition_Y']
            optDF = pd.DataFrame() 
            if col in nonIntCols:
                optDF[col] = df[col]
            else:
                df[col] = df[col].apply(np.log2)

        df = df.replace('', np.nan)
        df = df.replace(0, np.nan)  #Could this have any negative consequences?
        df = df.replace(float('-inf'), np.nan)
        df = df.replace(np.inf, np.nan)
        df = df.dropna() #default options should be appropriate

        #output normalized data:
        #If normalized file for well, open and append, otherwise new 
        well = fn.split('/')[-1].split('_')[0]
        #for fn2 in sorted(glob.glob(normalizedOutputFolder)):
            #if well in fn2:
            #    wellData=pd.read_table(fn)
            #    wellData.concat(df)
            #    wellData.to_csv(fn2)
            #else:
        filename = fn
        filename_out = normalizedOutputFolder + 'normalized' + filename.split('/')[-1]
        print(filename_out)
        df.to_csv(filename_out)

    return finaldf

#markerList = ['Rb-MAVS','LAMP2','M-IRF3','p-TBK1','NFAT-C1','NFkB','TMEM137','IRF5','COX4','STAT6','IRF7','IRF1','STAT5b','LC3A/B','STAT5a','Stat3','PKM2','p-mTOR']#,'p-p38 MAPK1','Stat1','p-Stat6']
##Is this right order for labeling below?
##Cycle 1: green (FITC?) - yellow (Cy3?) - red (Cy5?), Cycle 2 green - yellow - red, etc, doesn't specify 2o ab. 
#fn = '_B02_fld1_crop_cytoMasked.txt'
#numCycles=6

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







