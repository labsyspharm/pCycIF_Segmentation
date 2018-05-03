import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
from __future__ import division 

liveCellCount = []
deadCellCount = []
sorafLiveCellCount = []
sorafDeadCellCount = []
dmsoLiveCellCount = []
dmsoDeadCellCount = []



for fn in sorted(glob.glob('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/rawCSVs/*.csv')):
#    print(fn)
    df = pd.read_csv(fn)    
    if 'sig' in fn:
        nuclei = df.filter(items=['DNA-1_Nuc', 'DNA-2_Nuc', 'DNA-3_Nuc', 'DNA-4_Nuc','DNA-5_Nuc'])
    elif 'prot' in fn:
        nuclei = df.filter(items=['DNA-1_Nuc', 'DNA-2_Nuc', 'DNA-3_Nuc', 'DNA-4_Nuc'])
    else:
        print('fail')
        break
    model = KMeans(n_clusters=2,init='k-means++')
    model.fit(nuclei.values)  
    labels = model.labels_  
    clusters = pd.DataFrame(data=labels,columns=['cluster'])
    centers = model.cluster_centers_
    print(centers)
    liveCluster = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
    deadCluster = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])
    print(liveCluster)
    idx = clusters.index[clusters['cluster']==liveCluster].tolist()
    livecells = df.loc[idx]
    idx = clusters.index[clusters['cluster']==deadCluster].tolist()
    deadcells = df.loc[idx]

    liveCellCount.append(len(livecells))
    deadCellCount.append(len(deadcells))

    #split on drugs
    sorafLiveIdx = livecells.index[livecells['Drug']==7].tolist()
    sorafDiveIdx = deadcells.index[deadcells['Drug']==7].tolist()
    sorafLive = livecells.loc[sorafLiveIdx]
    sorafDead = deadcells.loc[sorafDiveIdx]

    dmsoLiveIdx = livecells.index[livecells['Drug']==11].tolist()
    dmsoDiveIdx = deadcells.index[deadcells['Drug']==11].tolist()
    dmsoLive = livecells.loc[dmsoLiveIdx]
    dmsoDead = deadcells.loc[dmsoDiveIdx]

    sorafLiveCellCount.append(len(sorafLive))
    sorafDeadCellCount.append(len(sorafDead))
    dmsoLiveCellCount.append(len(dmsoLive))
    dmsoDeadCellCount.append(len(dmsoDead))


totalCount = [deadCellCount + liveCellCount for deadCellCount, liveCellCount in zip(deadCellCount, liveCellCount)]
percentTotal = [deadCellCount/totalCount for deadCellCount, totalCount in zip(deadCellCount, totalCount)]
totalSoraf = [sorafDeadCellCount + sorafLiveCellCount for sorafDeadCellCount, sorafLiveCellCount in zip(sorafDeadCellCount, sorafLiveCellCount)]
percentSoraf = [sorafDeadCellCount/totalSoraf for sorafDeadCellCount, totalSoraf in zip(sorafDeadCellCount, totalSoraf)]
totaltDMSO = [dmsoDeadCellCount + dmsoLiveCellCount for dmsoDeadCellCount, dmsoLiveCellCount in zip(dmsoDeadCellCount, dmsoLiveCellCount)]
percentDMSO = [dmsoDeadCellCount/totaltDMSO for dmsoDeadCellCount, totaltDMSO in zip(dmsoDeadCellCount, totaltDMSO)]

columns = ['Total Live Cells', 'Total Lost Cells', 'Percent Lost Cells', 'Sorafenib Live Cells', 'Sorafenib Lost Cells', 'Sorafenib Percent Lost Cells', 'DMSO Live Cells', 'DMSO Lost Cells', 'DMSO Percent Lost Cells']
cellCountDF = pd.DataFrame(columns = columns)

cellCountDF['Total Live Cells'] = liveCellCount
cellCountDF['Total Lost Cells'] = deadCellCount
cellCountDF['Percent Lost Cells'] = percentTotal
cellCountDF['Sorafenib Live Cells'] = sorafLiveCellCount
cellCountDF['Sorafenib Lost Cells'] = sorafDeadCellCount
cellCountDF['Sorafenib Percent Lost Cells'] = percentSoraf
cellCountDF['DMSO Live Cells'] = dmsoLiveCellCount
cellCountDF['DMSO Lost Cells'] = dmsoDeadCellCount
cellCountDF['DMSO Percent Lost Cells'] = percentDMSO

cellCountDF.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/cellCount.csv')

    #write files
#    soraf.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/clusteredData/soraf/soraf_%s_%s' % (fn.split('_')[3], fn.split('_')[4]))
#    dmso.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/clusteredData/dmso/dmso_%s_%s' % (fn.split('_')[3], fn.split('_')[4]))



