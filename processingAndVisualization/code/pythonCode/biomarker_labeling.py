import pandas as pd
from copy import deepcopy

#Need checks for length based on cycle number, check new col length vs old. 
    

def labelBiomarkers(fn,numCycles,markerList=None,bleached=1):
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
#    outputDF.to_csv(newFN)

    return outputDF



