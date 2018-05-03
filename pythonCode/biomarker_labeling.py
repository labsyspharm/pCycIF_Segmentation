#Things need to be separated by antibody set. 
#B&C, D&E, F&G for proteome plates (6-8)
#Seemingly should be 4 nuclei, 12 cyto
#have 11 nuclei, 11 cyto

#For proteome there should be 16 markers, nuclear and cyto
#For signaling, 20.
#dapi 2-4, FITC 1-4, Cy3 1-4, Cy5 1-4
#Nuclear then cytoplasm

testPath = 'plate6_tif';
nucleiStack = [1 4];
row = 'B':'C';
col = 11;


#Its going to be important to separate ab sets since labels will be different
#ab sets could also define nuclear stack if we ignore plate 5
#fields for same well should be combined
#can combine different wells with same ab_set, but need to be able to note drug, dose, time 
#output file name: _B11_fl2_cytoMasked 
#keeps row and col so we can get drug and dose, doesn't keep plate currently.
#Should probably process on a whole plate level, maybe with option to drop parts
#This works well with current setup and normalization approach
