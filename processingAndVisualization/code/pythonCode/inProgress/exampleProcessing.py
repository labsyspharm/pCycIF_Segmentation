from __future__ import division 
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os
from cycifProcessingTools import thresholdNuclei
from cycifProcessingTools import normalizeNC


#Want to be more careful about file input/handling
#need to figure out where can and need to run code from, and to make sure python methods remain callable
outputDir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/output/normalized/'
for fn in sorted(glob.glob(outputDir+'norm*.txt')):
    newfn = thresholdNuclei(fn)
    finalDF = normalizeNC(newfn,df=None)


#input > output > ouput/normalization > lots of files. 
#Now automatically adding filders to output/ for all processing steps




#if __name__ == '__main__':
#    genes = ['EGF', 'EGFR', 'ERBB2', 'GRB2', 'SOS1', 'HRAS', 'RAF1',
#            'MAP2K1', 'MAPK1']

#    with open('reading/model.pkl', 'rb') as f:
#        model = pickle.load(f)
#    stmts1 = []
#    for k, v in model.items():
#        stmts1 += v
#    model, my_stmts3, my_stmts4 = get_subnetwork(stmts1, genes)


