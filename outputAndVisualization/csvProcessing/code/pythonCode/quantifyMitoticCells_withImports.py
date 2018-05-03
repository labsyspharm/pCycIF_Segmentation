from __future__ import division 
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os

import cycifProcessingTools as cpt

def checkClusterOutliers(df,idxHi,idxLo,pRbHiCells, pRbLoCells):
   
    #Handle outliers
    #if one cluster ends up with a very small number of cells
    #remove these cells and re-cluster   
    while len(idxHi) < 5:
        df = df.drop(df.index[idxHi])
        df = df.reset_index()

        pRbHiCells, pRbLoCells = clusterRb(fn,drugName,df)

    while len(idxLo) < 10:
        df = df.drop(df.index[idxHi])
        df = df.reset_index()

        pRbHiCells, pRbLoCells = clusterRb(fn,drugName,df)

    return pRbHiCells, pRbLoCells

def clusterRb(fn,df):

    #not sure if I still need this, but don't think it does any harm
    #Reset dataframe index after removing/slicing original df
    df = df.reset_index()

    if 'p-Rb_Nuc' in df.columns:
        pRb_Nuc = df['p-Rb_Nuc']
        pRb_Cyt = df['p-Rb_Cyto']
        pRb_NC = pRb_Nuc/pRb_Cyt

    elif 'pRb_Nuc' in df.columns:
        pRb_Nuc = df['pRb_Nuc']
        pRb_Cyt = df['pRb_Cyto']
        pRb_NC = pRb_Nuc/pRb_Cyt
    else:
        print('No pRb found in %s' % fn)

    #Try/Except handles files where there is no pRb
    #Enters NaN values 
    try:

        print(pRb_NC.shape)

        if pRb_NC.shape[0] > 1:
            model = KMeans(n_clusters=2,init='k-means++')
            model.fit(pRb_NC.values.reshape(-1,1))
            labels = model.labels_  
            clusters = pd.DataFrame(data=labels,columns=['cluster'])
            centers = model.cluster_centers_

            #print(centers)
            pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
            pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

            idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
            pRbHiCells = df.loc[idxHi]

            idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
            pRbLoCells = df.loc[idxLo]
        
            #Check for outliers, reclustering if needed
            pRbHiCells, pRbLoCells = checkClusterOutliers(df,idxHi,idxLo,pRbHiCells,pRbLoCells)
    
        else:
            pRbCount = [np.NaN, np.NaN]
            pRbHiCells = None
            pRbLoCells = None

    except NameError as err:
#        print(fn)
        pRbCount = [np.NaN, np.NaN]
        pRbHiCells = None
        pRbLoCells = None

    return pRbHiCells, pRbLoCells

def writeFiles(dfHi,dfLo,fn):

    #write files
    abset = fn.split('.')[0][:-1].split('_')[-1]
    drugDict = {2:'Afatinib',3:'Dasatinib',4:'Gefitinib',5:'Lapatinib',6:'Nilotinib',7:'Sorafenib',8:'Sunitinib',9:'Crizotinib',10:'Ponatinib',11:'DMSO'}
    drugName = drugDict[list(dfHi['Drug'])[0]]
    dirName = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/spreadsheets/highdose/%s/' % (drugName)
    if not os.path.exists(dirName):
        os.makedirs(dirName)    

    dirNameHi = dirName+'pRb_high/'
    dirNameLo = dirName+'pRb_low/'
    if not os.path.exists(dirNameHi):
        os.makedirs(dirNameHi)
    if not os.path.exists(dirNameLo):
        os.makedirs(dirNameLo)

    plate = fn.split('ProcessedData')[1][:2]

    dfHi.to_csv(dirNameHi+'%s_highdose_pRbHi_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))
    dfLo.to_csv(dirNameLo+'%s_lowdose_pRbLo_%s_%s.csv' % (drugName,fn.split('/')[-1].split('_')[-1].split('.')[0],plate))

    #count mitotic cells
    pRbCount = [dfHi.shape[0], dfLo.shape[0]]
    try:
        rowname = drugName+'_'+fn.split('_')[-1].split('.')[0]
        mitoticCellsDF.loc[rowname] = pRbCount
    except NameError:
        columns = ['pRb N/C High', 'pRb N/C Low']
        mitoticCellsDF = pd.DataFrame(columns = columns)
#        print(rowname)
#        print(mitoticCellsDF)
        rowname = drugName+'_'+fn.split('_')[-1].split('.')[0]

    mitoticCellsDF.loc[rowname] = pRbCount
    fn = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/pRbCellCount.csv'
#    with open(fn,'a') as f:
#        mitoticCellsDF.to_csv(f) 
    mitoticCellsDF.to_csv(fn, mode='a', header=False) 


foldername = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/rawCSVs/'
for fn in sorted(glob.glob('%s*.csv' % foldername)):
    if 'P3' not in fn:
#    if 'signaling2' in fn:
        print(fn)
        liveCellsDf = cpt.clusterNuclei(fn)
#        print('num live cells is %s' % liveCellsDf.shape[0])
        NC_DF = cpt.normalizeNC(fn,liveCellsDf)
        highDoseDF = cpt.splitDoses(fn,NC_DF)
        print('Num of high dose live cells for all drugs is %s' % highDoseDF.shape[0])
        drugDfList = cpt.splitDrugs(fn,highDoseDF,writeFile=1)
        i = 1
        for df in drugDfList:
            print(i)
            i=i+1
            print('Num of live cells for drug is %s' % df.shape[0])
            pRbHiCells, pRbLoCells = clusterRb(fn,df)
            if pRbHiCells is not None:
                writeFiles(pRbHiCells,pRbLoCells,fn)




#E,F,F 6





