from __future__ import division 
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import glob 
import os
#from cycifProcessingTools import thresholdNuclei
#from cycifProcessingTools import normalizeNC


#Want to be more careful about file input/handling
#need to figure out where can and need to run code from, and to make sure python methods remain callable
for fn in sorted(glob.glob('norm*.txt')):
    newfn = thresholdNuclei(fn)
    finalDF = normalizeNC(newfn,df=None)



