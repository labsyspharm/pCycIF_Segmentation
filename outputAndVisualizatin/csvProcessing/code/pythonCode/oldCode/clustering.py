import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 

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
    print(liveCluster)
    idx = clusters.index[clusters['cluster']==liveCluster].tolist()
    livecells = df.loc[idx]

    #split on drugs
    soraf_idx = livecells.index[livecells['Drug']==7].tolist()
    soraf = livecells.loc[soraf_idx]

    dmso_idx = livecells.index[livecells['Drug']==11].tolist()
    dmso = livecells.loc[dmso_idx]

    #write files
    soraf.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/clusteredData/soraf/soraf_%s_%s' % (fn.split('_')[3], fn.split('_')[4]))
    dmso.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/clusteredData/dmso/dmso_%s_%s' % (fn.split('_')[3], fn.split('_')[4]))



