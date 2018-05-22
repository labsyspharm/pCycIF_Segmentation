from __future__ import division 
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os


#def clusterNuclei(fn):

#    df = pd.read_csv(fn)    
#    colNames = list(df.columns)
#    nucleiCols = [col for col in colNames if 'DNA' in col and 'Nuc' in col]
#    nuclei = df.filter(items=nucleiCols)

#    model = KMeans(n_clusters=2,init='k-means++')
#    model.fit(nuclei.values)  
#    labels = model.labels_  
#    clusters = pd.DataFrame(data=labels,columns=['cluster'])
#    centers = model.cluster_centers_
#    liveCluster = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
#    idx = clusters.index[clusters['cluster']==liveCluster].tolist()
#    livecells = df.loc[idx]
#    
#    print('Num of live cells is %s' % livecells.shape[0])
#    print('Num of all cells is %s' % nuclei.shape[0])

#    #Should probably be able to output a file here
#    return livecells


def thresholdNuclei(fn):

    df = pd.read_csv(fn)    
    colNames = list(df.columns)
    nucleiCols = [col for col in colNames if 'DNA' in col and 'Nuc' in col]
    nuclei = df.filter(items=nucleiCols)
    keyVal = df[nucleiCols[0]].mean()*.5
    for col in nucleiCols:
        liveCells = df[df[col]>=keyVal]
        deadCells = df[df[col]<keyVal]

    print('Num of live cells is %s' % liveCells.shape[0])
    print('Num of lost cells is %s' % deadCells.shape[0])

    print(fn)
#    baseFolder = 
#    newFolder = #Add a ../ here? Eventually keep all output folders on same level, instead of a deep tree. Need to think about naming - use all transformations or only most recent?
#    if not os.path.exists(newFolder):
#        os.makedirs(newFolder)
#    outputFN = newFolder+fn
#    liveCells.to_csv(outputFN)

#    return outputFN

#Could also attempt boolean mask, maybe in numpy
#Write function and use  apply on all cols



def normalizeNC(fn,df=None):
    #Is it necessary to require df input? 
    nucColumns = []
    if not df:
        df = pd.read_csv(fn)   
    for col in df.columns:
        if '_Nuc' in col:
            nucColumns.append(col)

    for col in nucColumns:
        col_nuc = df[col]
        col_cyt_name = str(col.split('_Nuc')[:-1][0]) + '_Cyto'
        col_cyt = df[col_cyt_name]
        col_nc_ratio = col_nuc/col_cyt
        newname = str(col.split('_')[:-1][0]) + '_NC_Ratio'
        df[newname] = col_nc_ratio
#    if writeFile:
    df.to_csv('NCRatio_'+fn)
    return df



#More compelex file writing?

#        dirName = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/%s/highdose/%s/' % (abset, drug)
#        if not os.path.exists(dirName):
#            os.makedirs(dirName)
#        highdose.to_csv(dirName+'highdose_%s' % (fn.split('/')[-1]))






######################
#Less generalizable under here

#def splitDoses(fn,df):
#    if 'signaling' in fn:
#        highdose = 10
#    else:
#        if 'P6' in fn:
#            if 'proteome3' in fn:
#                highdose = 3
#            else:
#                highdose = 10
#        else:
#            highdose = 3

#    highdose_idx = df.index[df['Dose']==highdose].tolist()
#    highdose = df.loc[highdose_idx]

#        lowdose_idx = df.index[df['Dose']==1].tolist()
#        lowdose = df.loc[lowdose_idx]
#    
#    return highdose


#def splitDrugs(fn,df,writeFile=None):
#    drugIndexList = list(range(2,12))
#    drugDfList = []
#    drugDict = {2:'Afatinib',3:'Dasatinib',4:'Gefitinib',5:'Lapatinib',6:'Nilotinib',7:'Sorafenib',8:'Sunitinib',9:'Crizotinib',10:'Ponatinib',11:'DMSO'}
#    for drugID in drugIndexList:
#        drugName = drugDict[drugID]
#        drug_idx = df.index[df['Drug']==drugID].tolist()   
#        drug_df = df.loc[drug_idx]

#        dose_df = splitDoses(fn,drugName,drug_df)
#        drugDfList.append(drug_df)

#        if writeFile:
#            abset = fn.split('.')[0][:-1].split('_')[-1]
#            dirName = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/%s/highdose/%s/' % (abset, drugName)
#            print(dirName)
#            if not os.path.exists(dirName):
#                os.makedirs(dirName)
#            drug_df.to_csv(dirName+'%s_highdose_%s' % (drugName, fn.split('/')[-1]) )
#    return drugDfList



#def clusterObservable(fn,df,observable):

#    #not sure if I still need this, but don't think it does any harm
#    #Reset dataframe index after removing/slicing original df
#    df = df.reset_index()

#    if 'p-Rb_Nuc' in df.columns:
#        pRb_Nuc = df['p-Rb_Nuc']
#        pRb_Cyt = df['p-Rb_Cyto']
#        pRb_NC = pRb_Nuc/pRb_Cyt

#    elif 'pRb_Nuc' in df.columns:
#        pRb_Nuc = df['pRb_Nuc']
#        pRb_Cyt = df['pRb_Cyto']
#        pRb_NC = pRb_Nuc/pRb_Cyt
#    else:
#        print('No pRb found in %s' % fn)

#    #Try/Except handles files where there is no pRb
#    #Enters NaN values 
#    try:

#        model = KMeans(n_clusters=2,init='k-means++')
#        model.fit(pRb_Nuc.values.reshape(-1,1))
#        labels = model.labels_  
#        clusters = pd.DataFrame(data=labels,columns=['cluster'])
#        centers = model.cluster_centers_

#        #print(centers)
#        pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
#        pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

#        idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
#        pRbHiCells = df.loc[idxHi]

#        idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
#        pRbLoCells = df.loc[idxLo]
#    
#        #Check for outliers, reclustering if needed
#        pRbHiCells, pRbLoCells = checkClusterOutliers(df,idxHi,idxLo,pRbHiCells,pRbLoCells)

#    except NameError as err:
##        print(fn)
#        pRbCount = [np.NaN, np.NaN]
#        pRbHiCells = None
#        pRbLoCells = None

#    return pRbHiCells, pRbLoCells
