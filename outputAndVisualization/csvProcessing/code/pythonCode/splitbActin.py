from __future__ import division
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os

import cycifProcessingTools as cpt

def checkClusterOutliers(df,idxHi,idxLo,bActinHiCells, bActinLoCells):
   
    #Handle outliers
    #if one cluster ends up with a very small number of cells
    #remove these cells and re-cluster   
#    print('High Cell Count is %s' % len(idxHi))
#    print('Low Cell Count is %s' % len(idxLo))
    if len(idxHi) < 10:
        df = df.drop(df.index[idxHi])
        df = df.reset_index(drop=True)

        bActinHiCells, bActinLoCells = clusterbActin(fn,df)

    if len(idxLo) < 10:
        df = df.drop(df.index[idxLo])
        df = df.reset_index(drop=True)

        bActinHiCells, bActinLoCells = clusterbActin(fn,df)

    return bActinHiCells, bActinLoCells

def clusterbActin(fn,df):

    #not sure if I still need this, but don't think it does any harm
    #Reset dataframe index after removing/slicing original df
    df = df.reset_index(drop=True)

    if 'B-actin_Nuc' in df.columns:
        bActin_Nuc = df['B-actin_Nuc']
        bActin_Cyto = df['B-actin_Cyto']
        bActin_NC = bActin_Nuc/bActin_Cyto
    else:
        print('No bActin found in %s' % fn)

    #Try/Except handles files where there is no pRb
    #Enters NaN values 
    try:

        model = KMeans(n_clusters=2,init='k-means++')
        model.fit(bActin_Cyto.values.reshape(-1,1))
        labels = model.labels_  
        clusters = pd.DataFrame(data=labels,columns=['cluster'])
        centers = model.cluster_centers_

        #print(centers)
        bActinHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
        bActinLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

        idxHi = clusters.index[clusters['cluster']==bActinHi].tolist()
        bActinHiCells = df.loc[idxHi]

        idxLo = clusters.index[clusters['cluster']==bActinLo].tolist()
        bActinLoCells = df.loc[idxLo]
    
        #Check for outliers, reclustering if needed
        bActinHiCells, bActinLoCells = checkClusterOutliers(df,idxHi,idxLo,bActinHiCells,bActinLoCells)

    except NameError as err:
#        print(fn)
        #pRbCount = [np.NaN, np.NaN]
        bActinHiCells = None
        bActinLoCells = None

    return bActinHiCells, bActinLoCells

def writeFiles(dfHi,dfLo,fn):

    #write files
    abset = fn.split('.')[0][:-1].split('_')[-1]
    drugDict = {2:'Afatinib',3:'Dasatinib',4:'Gefitinib',5:'Lapatinib',6:'Nilotinib',7:'Sorafenib',8:'Sunitinib',9:'Crizotinib',10:'Ponatinib',11:'DMSO'}
    drugName = drugDict[list(dfHi['Drug'])[0]]
    dirName = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/bActinSplit/spreadsheets/highdose/%s/' % (drugName)
    if not os.path.exists(dirName):
        os.makedirs(dirName)    

    dirNameHi = dirName+'bActin_high/'
    dirNameLo = dirName+'bActin_low/'
    if not os.path.exists(dirNameHi):
        os.makedirs(dirNameHi)
    if not os.path.exists(dirNameLo):
        os.makedirs(dirNameLo)

    plate = fn.split('ProcessedData')[1][:2]
    dfHi.to_csv(dirNameHi+'%s_highdose_bActinHi_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))
    dfLo.to_csv(dirNameLo+'%s_lowdose_bActinLo_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))



foldername = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/rawCSVs/'
for fn in sorted(glob.glob('%s/*.csv' % foldername)):
    if 'old' not in fn:
        liveCellsDf = cpt.clusterNuclei(fn)
        NC_DF = cpt.normalizeNC(fn,liveCellsDf)
        highDoseDF = cpt.splitDoses(fn,NC_DF)
        drugDfList = cpt.splitDrugs(fn,highDoseDF)

        for df in drugDfList:
            bActinHiCells, bActinLoCells = clusterbActin(fn,df)
            if bActinHiCells is not None:
                writeFiles(bActinHiCells,bActinLoCells,fn)
















