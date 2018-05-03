#!/usr/bin/python
#TODO:
#Allow for more general channel/biomarker declaration
#pick dynamic max value for re-normalization

#Check matlab/python versions and dependencies on o2 - try version with no matlab binary

#WORKING MATLAB CALL
#./run_segmentScript.sh /o2/path/MATLAB_Runtime/v901/ '/image/files/location' '[1 cycleNumber]' 'rowLetters' '[columnNumbers]'

#WORKING PYTHON CALL
#subprocess.call(["./run_bobbySegmentCardiomyo.sh", "/home/bobby/MATLAB_Runtime/v901/", '/home/bobby/Dropbox/MATLAB/cardiotox_cycif/segmentation/bobbySegmentCardioMyo/plate6_tif', '[1 4]', "'BEG'", '[10:11]'])

import pandas as pd
import argparse
import glob
import sys 
import numpy as np
import subprocess
from ab_dict_signaling_file import ab_dict_signaling
from ab_dict_proteome_file import ab_dict_proteome


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


def resegement_plate(plate):
    platename = plate
    platenumber = plate[-1]
    if platenumber in ['1', '2', '3', '4', '5']:
        nuclearstack = '[1 5]'
    else:
        nuclearstack = '[1 4]'
    rows = 'BCDEFG'
    cols = '[2:11]'
    subprocess.call(["./run_bobbySegmentCardiomyo_orchestra.sh", "/home/rps21/MCR/v901/", '/n/scratch2/rps21/ResegmentedPlates/%s' % platename, nuclearstack, rows, cols])
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
            df = df.dropna() #default options should be appropriate

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
                    df[col] = df[col].apply(np.log2)

            df = df.replace('', np.nan)
            df = df.replace(0, np.nan)  #Could this have any negative consequences?
            df = df.replace(float('-inf'), np.nan)
            df = df.replace(np.inf, np.nan)
            df = df.dropna() #default options should be appropriate


            #can i automate building this list based on number of cycles?
            proteome_channels = ['DAPI-0001','DAPI-0002','DAPI-0003','DAPI-0004','FITC-0001','FITC-0002','FITC-0003','FITC-0004','Cy3-0001','Cy3-0002','Cy3-0003','Cy3-0004','Cy5-0001','Cy5-0002','Cy5-0003','Cy5-0004']
            signaling_channels = ['DAPI-0001','DAPI-0002','DAPI-0003','DAPI-0004','DAPI-0005','FITC-0001','FITC-0002','FITC-0003','FITC-0004','FITC-0005','Cy3-0001','Cy3-0002','Cy3-0003','Cy3-0004','Cy3-0005','Cy5-0001','Cy5-0002','Cy5-0003','Cy5-0004','Cy5-0005']


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

            new_cols = df.columns + optDF.columns
            df.columns = new_cols
 
            #add dose and drug columns
            doseDrug = filename.split('_')[1]
            row = doseDrug[0] #Need to convert string to num into final csv
            if plate == 'ProcessedDataP6':
                dose_dict = {('B','proteome1'):10,('D','proteome2'):1,('F','proteome3'):3,('C','proteome1'):3,('E','proteome2'):10,('G','proteome3'):1}
            else:
                dose_dict = {('B','signaling1'):10,('E','signaling2'):10,('C','signaling1'):3,('F','signaling2'):3,('B','proteome1'):3,('D','proteome2'):3,('F','proteome3'):3,('D','signaling1'):1,('G','signaling2'):1,('C','proteome1'):1,('E','proteome2'):1,('G','proteome3'):1}
            dose = dose_dict[(row,ab_set)]
            df['Dose'] = pd.Series([dose]*len(df.index),index = df.index)
    
            drug = doseDrug[1:]
            df['Drug'] = pd.Series([drug]*len(df.index),index = df.index)

            if filenum == 1:
                filenum = 2
                finaldf = df
            else:
                finaldf = pd.concat([finaldf,df],axis=0,ignore_index=True)

#    filename_out = 'combined_abset.csv' #THIS IS GOING TO OVERWRITE
    filename_out = 'resegmented_normalized_logtransformed_%s_%s.csv' % (plate, ab_set) 
    finaldf.to_csv(filename_out)


#plates_to_use = ['plate1','plate6'] #Hard code for now
##plates_to_use = ['plate6'] #Hard code for now
#for plate in plates_to_use:
#resegement_plate(plate)

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


