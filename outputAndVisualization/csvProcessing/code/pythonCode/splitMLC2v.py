from __future__ import division 
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os

import cycifProcessingTools as cpt

def checkClusterOutliers(df,idxHi,idxLo,mlc2vHiCells, mlc2vLoCells):
   
    #Handle outliers
    #if one cluster ends up with a very small number of cells
    #remove these cells and re-cluster   
    while len(idxHi) < 10:
        df = df.drop(df.index[idxHi])
        df = df.reset_index()

        mlc2vHiCells, mlc2vLoCells = clusterMlc2v(fn,drugName,df)

    while len(idxLo) < 10:
        df = df.drop(df.index[idxHi])
        df = df.reset_index()

        mlc2vHiCells, mlc2vLoCells = clusterMlc2v(fn,df)

    return mlc2vHiCells, mlc2vLoCells

def clusterMlc2v(fn,df):

    #not sure if I still need this, but don't think it does any harm
    #Reset dataframe index after removing/slicing original df
    df = df.reset_index()

    if 'MLC2v_Nuc' in df.columns:
        mlc2v_Nuc = df['MLC2v_Nuc']
        mlc2v_Cyto = df['MLC2v_Cyto']
        mlc2v_NC = mlc2v_Nuc/mlc2v_Cyto
    else:
        print('No MLC2v found in %s' % fn)

    #Try/Except handles files where there is no pRb
    #Enters NaN values 
    try:

        model = KMeans(n_clusters=2,init='k-means++')
        model.fit(mlc2v_Cyto.values.reshape(-1,1))
        labels = model.labels_  
        clusters = pd.DataFrame(data=labels,columns=['cluster'])
        centers = model.cluster_centers_

        #print(centers)
        mlc2vHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
        mlc2vLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

        idxHi = clusters.index[clusters['cluster']==mlc2vHi].tolist()
        mlc2vHiCells = df.loc[idxHi]

        idxLo = clusters.index[clusters['cluster']==mlc2vLo].tolist()
        mlc2vLoCells = df.loc[idxLo]
    
        #Check for outliers, reclustering if needed
        mlc2vHiCells, mlc2vLoCells = checkClusterOutliers(df,idxHi,idxLo,mlc2vHiCells,mlc2vLoCells)

    except NameError as err:
#        print(fn)
        #pRbCount = [np.NaN, np.NaN]
        mlc2vHiCells = None
        mlc2vLoCells = None

    return mlc2vHiCells, mlc2vLoCells

def writeFiles(dfHi,dfLo,fn):

    #write files
    abset = fn.split('.')[0][:-1].split('_')[-1]
    drugDict = {2:'Afatinib',3:'Dasatinib',4:'Gefitinib',5:'Lapatinib',6:'Nilotinib',7:'Sorafenib',8:'Sunitinib',9:'Crizotinib',10:'Ponatinib',11:'DMSO'}
    drugName = drugDict[list(dfHi['Drug'])[0]]
    dirName = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/mlc2vSplit/spreadsheets/highdose/%s/' % (drugName)
    if not os.path.exists(dirName):
        os.makedirs(dirName)    

    dirNameHi = dirName+'mlc2v_high/'
    dirNameLo = dirName+'mlc2v_low/'
    if not os.path.exists(dirNameHi):
        os.makedirs(dirNameHi)
    if not os.path.exists(dirNameLo):
        os.makedirs(dirNameLo)

    plate = fn.split('ProcessedData')[1][:2]
    dfHi.to_csv(dirNameHi+'%s_highdose_mlc2vHi_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))
    dfLo.to_csv(dirNameLo+'%s_lowdose_mlc2vLo_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))



foldername = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/rawCSVs/'
for fn in sorted(glob.glob('%s/*.csv' % foldername)):
    if 'old' not in fn:
        liveCellsDf = cpt.clusterNuclei(fn)
        NC_DF = cpt.normalizeNC(fn,liveCellsDf)
        highDoseDF = cpt.splitDoses(fn,NC_DF)
        drugDfList = cpt.splitDrugs(fn,highDoseDF)

        for df in drugDfList:
            mlc2vHiCells, mlc2vLoCells = clusterMlc2v(fn,df)
            if mlc2vHiCells is not None:
                writeFiles(mlc2vHiCells,mlc2vLoCells,fn)











