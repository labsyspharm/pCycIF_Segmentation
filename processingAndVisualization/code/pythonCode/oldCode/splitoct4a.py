from __future__ import division 
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os

import cycifProcessingTools as cpt

def checkClusterOutliers(df,idxHi,idxLo,oct4aHiCells, oct4aLoCells):
   
    #Handle outliers
    #if one cluster ends up with a very small number of cells
    #remove these cells and re-cluster   
#    print('High Cell Count is %s' % len(idxHi))
#    print('Low Cell Count is %s' % len(idxLo))
    if len(idxHi) < 10:
        df = df.drop(df.index[idxHi])
        df = df.reset_index(drop=True)

        oct4aHiCells, oct4aLoCells = clusterOct4a(fn,df)

    if len(idxLo) < 10:
        df = df.drop(df.index[idxLo])
        df = df.reset_index(drop=True)

        oct4aHiCells, oct4aLoCells = clusterOct4a(fn,df)

    return oct4aHiCells, oct4aLoCells

def clusterOct4a(fn,df):

    #not sure if I still need this, but don't think it does any harm
    #Reset dataframe index after removing/slicing original df
    df = df.reset_index(drop=True)

    if 'Oct4a_Nuc' in df.columns:
        oct4a_Nuc = df['Oct4a_Nuc']
        oct4a_Cyto = df['Oct4a_Cyto']
        oct4a_NC = oct4a_Nuc/oct4a_Cyto
    else:
        print('No oct4a found in %s' % fn)

    #Try/Except handles files where there is no pRb
    #Enters NaN values 
    try:

        model = KMeans(n_clusters=2,init='k-means++')
        model.fit(oct4a_Cyto.values.reshape(-1,1))
        labels = model.labels_  
        clusters = pd.DataFrame(data=labels,columns=['cluster'])
        centers = model.cluster_centers_

        #print(centers)
        oct4aHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
        oct4aLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

        idxHi = clusters.index[clusters['cluster']==oct4aHi].tolist()
        oct4aHiCells = df.loc[idxHi]

        idxLo = clusters.index[clusters['cluster']==oct4aLo].tolist()
        oct4aLoCells = df.loc[idxLo]
    
        #Check for outliers, reclustering if needed
        oct4aHiCells, oct4aLoCells = checkClusterOutliers(df,idxHi,idxLo,oct4aHiCells,oct4aLoCells)

    except NameError as err:
#        print(fn)
        #pRbCount = [np.NaN, np.NaN]
        oct4aHiCells = None
        oct4aLoCells = None

    return oct4aHiCells, oct4aLoCells

def writeFiles(dfHi,dfLo,fn):

    #write files
    abset = fn.split('.')[0][:-1].split('_')[-1]
    drugDict = {2:'Afatinib',3:'Dasatinib',4:'Gefitinib',5:'Lapatinib',6:'Nilotinib',7:'Sorafenib',8:'Sunitinib',9:'Crizotinib',10:'Ponatinib',11:'DMSO'}
    drugName = drugDict[list(dfHi['Drug'])[0]]
    dirName = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/oct4aSplit/spreadsheets/highdose/%s/' % (drugName)
    if not os.path.exists(dirName):
        os.makedirs(dirName)    

    dirNameHi = dirName+'oct4a_high/'
    dirNameLo = dirName+'oct4a_low/'
    if not os.path.exists(dirNameHi):
        os.makedirs(dirNameHi)
    if not os.path.exists(dirNameLo):
        os.makedirs(dirNameLo)

    plate = fn.split('ProcessedData')[1][:2]
    dfHi.to_csv(dirNameHi+'%s_highdose_oct4aHi_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))
    dfLo.to_csv(dirNameLo+'%s_lowdose_oct4aLo_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))



foldername = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/rawCSVs/'
for fn in sorted(glob.glob('%s/*.csv' % foldername)):
    if 'old' not in fn:
        liveCellsDf = cpt.clusterNuclei(fn)
        NC_DF = cpt.normalizeNC(fn,liveCellsDf)
        highDoseDF = cpt.splitDoses(fn,NC_DF)
        drugDfList = cpt.splitDrugs(fn,highDoseDF)

        for df in drugDfList:
            oct4aHiCells, oct4aLoCells = clusterOct4a(fn,df)
            if oct4aHiCells is not None:
                writeFiles(oct4aHiCells,oct4aLoCells,fn)











