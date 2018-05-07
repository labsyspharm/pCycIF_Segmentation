import pandas as pd
from copy import deepcopy

#Need checks for length based on cycle number, check new col length vs old. 
    

def labelBiomarkers(fn,markerList=None,bleached=1):
    columnLabels = []
    numMarkers = numCycles*4
    dapiList = []


    #If no list of markers is provided, create generic naming
    if markerList:
        #add dapi labels to marker list. May need to handle this being optional 
        for i in range(0,numCycles):
            dapiList.append('DNA-%s' % str(i+1))
        markerList = dapiList+markerList
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_Nuc')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_Nuc')
        for i in range(0,len(markerList)):
            columnLabels.append(markerList[i]+'_Cyto')
            if bleached == 1 and (i+1) % numCycles != 0:   
                columnLabels.append('Bleached-%s' % markerList[i] + '_Cyto')


    #If list is provided, use for labeling
    else:
        channelList = ['DAPI','FITC','Cy3','Cy5']
        for i in range(1,numCycles+1):
            for i in range(0,numCycles):
                dapiList.append('DNA-%s' % str(i+1) + '_Nuc')
            for channel in channelList:
                label ='-%s' % i
                columnLabels.append(label + '_Nuc')
        for i in range(1,numCycles+1):
            for i in range(0,numCycles):
                dapiList.append('DNA-%s' % str(i+1) + '_Cyto')
            for channel in channelList:
                label ='-%s' % i
                columnLabels.append(label + '_Cyto')
    #missing and bleach 

    #For now hard code, may want to read these from original file, or get input from user
    extraColumns = ['NucleusArea','CytoplasmArea','CellPosition_X','CellPosition_Y']
    columnLabels = columnLabels + extraColumns

    originalDF = pd.read_table(fn)
    outputDF = deepcopy(originalDF)
    outputDF.columns = columnLabels
    newFN = fn.split('.')[0]+'.csv'
    outputDF.to_csv(newFN)

    return outputDF


##########
##Testing
#markerList = ['Rb-MAVS','LAMP2','M-IRF3','p-TBK1','NFAT-C1','NFkB','TMEM137','IRF5','COX4','STAT6','IRF7','IRF1','STAT5b','LC3A/B','STAT5a','Stat3','PKM2','p-mTOR']#,'p-p38 MAPK1','Stat1','p-Stat6']
##Is this right order for labeling below?
##Cycle 1: green (FITC?) - yellow (Cy3?) - red (Cy5?), Cycle 2 green - yellow - red, etc, doesn't specify 2o ab. 
#fn = '_B02_fld1_crop_cytoMasked.txt'
#numCycles=6

#df = labelBiomarkers(fn,markerList,1)

