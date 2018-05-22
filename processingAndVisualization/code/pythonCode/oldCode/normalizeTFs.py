import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 

TFlist = ['p-c-Jun_Cyto','p-STAT5_Cyto','p-STAT1_Cyto','p-STAT3_Cyto','c-Myc_Cyto','c-Jun_Cyto','NFkB_Cyto','p-STAT5_Cyto','STAT3_Cyto']
for drug in ['dmso','soraf']:
    for fn in sorted(glob.glob('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/proteome/**/%s/*.csv' % drug,recursive=True)):
        print(fn)
        df = pd.read_csv(fn)
        for tf in TFlist:
            if tf in df.columns:
                tf_cyto = df[tf]
                tf_nuc_name = str(tf.split('_')[:-1][0]) + '_Nuc'
                tf_nuc = df[tf_nuc_name]
                tf_nc_ratio = tf_nuc/tf_cyto
                newname = str(tf.split('_')[:-1][0]) + '_NC_Ratio'
                df[newname] = tf_nc_ratio
                #pull tf_cyto and nuc 
        dose = fn.split('/')[-1].split('_')[0]
        df.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/proteome/%s/%s/%s' % (dose,drug,fn.split('/')[-1]))



for drug in ['dmso','soraf']:
    for fn in sorted(glob.glob('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/doseSeparated/signaling/**/%s/*.csv' % drug,recursive=True)):
        df = pd.read_csv(fn)
        for tf in TFlist:
            if tf in df.columns:
                tf_cyto = df[tf]
                tf_nuc_name = str(tf.split('_')[:-1][0]) + '_Nuc'
                tf_nuc = df[tf_nuc_name]
                tf_nc_ratio = tf_nuc/tf_cyto
                newname = str(tf.split('_')[:-1][0]) + '_NC_Ratio'
                df[newname] = tf_nc_ratio
                #pull tf_cyto and nuc 
        dose = fn.split('/')[-1].split('_')[0]
        df.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/%s/%s/%s' % (dose,drug,fn.split('/')[-1]))







#        lowdose.to_csv('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/clustering/doses/signaling/lowdose/%s/lowdose_%s' % (drug,fn.split('/')[-1]))
