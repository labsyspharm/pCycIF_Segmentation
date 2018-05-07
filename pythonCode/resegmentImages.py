#!/usr/bin/python
#TODO:
#Allow for more general channel/biomarker declaration
#pick dynamic max value for re-normalization

#Check matlab dependencies on orchestra - try matlab binary without python first
#check python versions/dependencies

#WORKING MATLAB CALL
#./run_bobbySegmentCardiomyo.sh /home/bobby/MATLAB_Runtime/v901/ '/
#home/bobby/Dropbox/MATLAB/cardiotox_cycif/segmentation/bobbySegmentCardioMyo/plate6_tif' '[1 4]' 'BEG' '[10:11]'

#WORKING PYTHON CALL
#subprocess.call(["./home/rps21/Cardiomyocite/segmentation/bobbySegmentCardioMyo/run_bobbySegmentCardiomyo.sh", "/home/bobby/MATLAB_Runtime/v901/", '/home/bobby/Dropbox/MATLAB/cardiotox_cycif/segmentation/bobbySegmentCardioMyo/plate6_tif', '[1 4]', "'BEG'", '[10:11]'])



#python python_analysis_scripts/create_custom_csv.py --drugs Afatinib Sorafenib DMSO --dose 1 3.16 --time 1 4 24 48 --abs signaling1 signaling2 proteome1 proteome2 proteome3
#python python_analysis_scripts/create_custom_csv.py --drugs drug1 drug2 ... --dose dose1 dose2 ... --time time1 time2 ... --abs abset1 abset2 ...


import pandas as pd
import argparse
import glob
import sys 
from ab_dict_signaling_file import ab_dict_signaling
from ab_dict_proteome_file import ab_dict_proteome
#from drug_code_dict import drug_code
import numpy as np
import subprocess

#Take plate folder as input 
platechoices = ['ProcessedDataP1','ProcessedDataP2','ProcessedDataP3','ProcessedDataP4','ProcessedDataP5','ProcessedDataP6','ProcessedDataP7','ProcessedDataP8']
parser = argparse.ArgumentParser(description='Write')
parser.add_argument('--plate', type=str, nargs='+', help='<Required> Plates to resegment '+', '.join(platechoices), metavar='plate names', choices=platechoices, required=True)    
args = parser.parse_args()
platename = args.plate[0]


##Function to filter rows, plates, based on ab_set
##Call withing main function
def filter_rows(ab_set):

    if ab_set == 'signaling1':
        rowlist = ['B','C','D']
    elif ab_set == 'signaling2':
        rowlist = ['E','F','G']
    elif ab_set == 'proteome1':
        rowlist = ['B','C']
    elif ab_set == 'proteome2':
        rowlist = ['D','E']
    elif ab_set == 'proteome3':
        rowlist = ['F','G']
    rows_to_use = list(set(rowlist))
    return rows_to_use

#Assume we've run segmentation on a plate, have a folder of .txt files. Don't have to filter plate but will want some way to input/save info on plate/time

#This works, takes Clarence output, normalizes, combines plate 
#Want a wrapper to matlab to produce this output
#Potentially functionalize matlab script to allow input variables to define nuclear length and cols
#want to be able to script around that to cover all plates 
#Want to be able to further slice data so don't have include all col/rows/etc


#Matlab call
#subprocess.call(["./home/rps21/Cardiomyocite/segmentation/bobbySegmentCardioMyo/run_bobbySegmentCardiomyo.sh", "/home/bobby/MATLAB_Runtime/v901/", '/home/bobby/Dropbox/MATLAB/cardiotox_cycif/segmentation/bobbySegmentCardioMyo/%s' % platename, '[1 4]', "'BEG'", '[10:11]'])

#For available plates 
#

#platename = 'plate6'
#nuclearstack = '[1 4]'
#rows = 'BCDEFG'
#cols = '[10:11]'
#subprocess.call(["./run_bobbySegmentCardiomyo_orchestra.sh", "/home/rps21/MCR/v901/", '/home/rps21/Cardiomyocite/segmentation/bobbySegmentCardioMyo/%s' % platename, nuclearstack, rows, cols])
#This should build all the txt files we will use for normalization






def resegment_plate(plate):
    platename = plate
    platenumber = plate[-1]
    if platenumber in ['1', '2', '3', '4', '5']:
        nuclearstack = '[1 5]'
    else:
        nuclearstack = '[1 4]'
    print(platenumber)
    print(nuclearstack)
    rows = 'BCDEFG'
    cols = '[2:11]'
    print('/n/scratch2/rps21/ResegmentedPlates/%s' % platename, nuclearstack, rows, cols)
    subprocess.call(["/home/rps21/Cardiomyocite/segmentation/bobbySegmentCardioMyo/run_bobbySegmentCardiomyo_orchestra_addbackdapi.sh", "/home/rps21/MCR/v901/", '/n/scratch2/rps21/ResegmentedPlates/%s' % platename, nuclearstack, rows, cols])

#    print('done segmenting %s' % plate)


def normalize_plate(plate,ab_set):
    rows = filter_rows(ab_set)
    filenum = 1
    print('rows are %s' % rows)
    for row in rows:
         for fn in sorted(glob.glob('/n/scratch2/rps21/ResegmentedPlates/'+plate+'/_'+row+'*.txt')):
            filename = fn
            #for file, add column with dose and drug. Should have time but will be trickier
            #concat all files for this ab set and plate vertically
            df = pd.read_table(fn)
            #Normalize data
            for col in df.columns:
                newmin = min(df[col])
                newmax = 15000  #Really need a better way of setting this
                oldmax = max(df[col])
                if (oldmax - newmin) == 0:
                    print(fn)
                    df[col] = np.NaN
                else:
                    df[col] = (df[col] - newmin)*((newmax-newmin)/(oldmax - newmin))+newmin
                    df[col] = df[col].apply(np.log2)

            df = df.replace('', np.nan)
            df = df.replace(float('-inf'), np.nan)
            df = df.replace(np.inf, np.nan)
            df = df.dropna() #default options should be appropriate

            #can i automate building this list based on number of cycles?
            proteome_channels = ['FITC-0001','FITC-0002','FITC-0003','FITC-0004','Cy3-0001','Cy3-0002','Cy3-0003','Cy3-0004','Cy5-0001','Cy5-0002','Cy5-0003','Cy5-0004']
            signaling_channels = ['FITC-0001','FITC-0002','FITC-0003','FITC-0004','FITC-0005','Cy3-0001','Cy3-0002','Cy3-0003','Cy3-0004','Cy3-0005','Cy5-0001','Cy5-0002','Cy5-0003','Cy5-0004','Cy5-0005']

            new_cols = []
            if ab_set in ['signaling1', 'signaling2']:
                for channel in signaling_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_signaling[ab]+'_Nuc') #This file needs to be moved and imported
                for channel in signaling_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_signaling[ab]+'_Cyto') #This file needs to be moved and imported
            else:
                for channel in proteome_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_proteome[ab]+'_Nuc') #This file needs to be moved and imported
                for channel in proteome_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_proteome[ab]+'_Cyto') #This file needs to be moved and imported



            df.columns = new_cols

            #add dose and drug column s
            doseDrug = filename.split('_')[1]
            dose = doseDrug[0] #Need to convert string to go into final csv
            drug = doseDrug[1:]
            df['Drug'] = pd.Series([drug]*len(df.index),index = df.index)
            if filenum == 1:
                filenum = 2
                finaldf = df
            else:
                finaldf = pd.concat([finaldf,df],axis=0,ignore_index=True)

#    filename_out = 'combined_abset.csv' #THIS IS GOING TO OVERWRITE
    filename_out = 'resegmented_%s_%s.csv' % (plate, ab_set) 
    finaldf.to_csv(filename_out)


sigsets = ['signaling1','signaling2']
protsets = ['proteome1','proteome2','proteome3']
resegment_plate(platename)


#plates_to_use = ['plate1','plate6'] #Hard code for now
##plates_to_use = ['plate6'] #Hard code for now
#for plate in plates_to_use:
#    resegement_plate(plate)
#    if plate == 'plate1':
#        ab_sets = sigsets
#    else:
#        ab_sets = protsets
##    print(ab_sets)
#    for ab in ab_sets:
#        print('calling normalize with plate %s and ab_set %s' % (plate,ab))
#        normalize_plate(plate,ab)

#Normalize
sigsets = ['signaling1','signaling2']
protsets = ['proteome1','proteome2','proteome3']
platenumber = platename[-1]
if platenumber in ['1', '2', '3', '4', '5']:
    ab_sets = sigsets
else:
    ab_sets = protsets

for ab in ab_sets:
    #print('calling normalize with plate %s and ab_set %s' % (,ab))
    normalize_plate(platename,ab)

