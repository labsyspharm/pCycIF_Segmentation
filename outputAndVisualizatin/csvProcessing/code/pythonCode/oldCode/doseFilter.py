import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 

for drug in ['dmso','soraf']:
    for fn in sorted(glob.glob('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/clusteredData/proteome/%s/*.csv' % drug)):
    
#        print(fn)
        df = pd.read_csv(fn)    

        #split on drugs
        highdose_idx = df.index[df['Dose']==3].tolist()
        highdose = df.loc[highdose_idx]

        lowdose_idx = df.index[df['Dose']==1].tolist()
        lowdose = df.loc[lowdose_idx]

        #write files
        highdose.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/proteome/highdose/%s/highdose_%s' % (drug,fn.split('/')[-1]))
        lowdose.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/proteome/lowdose/%s/lowdose_%s' % (drug,fn.split('/')[-1]))


for drug in ['dmso','soraf']:
    for fn in sorted(glob.glob('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/clusteredData/signaling/%s/*.csv' % drug)):
#        print(fn)
        df = pd.read_csv(fn)    

        #split on drugs
        highdose_idx = df.index[df['Dose']==10].tolist()
        highdose = df.loc[highdose_idx]

        meddose_idx = df.index[df['Dose']==3].tolist()
        meddose = df.loc[highdose_idx]

        lowdose_idx = df.index[df['Dose']==1].tolist()
        lowdose = df.loc[lowdose_idx]

        #write files
        highdose.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/signaling/highdose/%s/highdose_%s' % (drug,fn.split('/')[-1]))
        meddose.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/signaling/meddose/%s/meddose_%s' % (drug,fn.split('/')[-1]))
        lowdose.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/signaling/lowdose/%s/lowdose_%s' % (drug,fn.split('/')[-1]))


