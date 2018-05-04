#5/4/18
#Can easily build the list/order of channels from number of cycles. 
#Need a list of markers with their channel/cycle 
#Could enforce an order, so that simply a correctly ordered list of biomarkers would work
#i.e. Cy3 1:end, Cy5 1:end, etc or Cy3 1, Cy5 1, FITC 1; Cy3 2, Cy 5 2, FITC 2; etc

#CHECK order in Clarence's output is still the same
#Is there a way to check for this? Check if first couple columns match expected
#If not, print warning and don't relabel 

#Techinically with the right order I can just loop through columns and replace
#But if I get them in an order I associate with channel/cycle, maybe I can handle unpredicted ordering from segmentation



 #can i automate building this list based on number of cycles?
            proteome_channels = ['DAPI-0001','DAPI-0002','DAPI-0003','DAPI-0004','FITC-0001','FITC-0002','FITC-0003','FITC-0004','Cy3-0001','Cy3-0002','Cy3-0003','Cy3-0004','Cy5-0001','Cy5-0002','Cy5-0003','Cy5-0004']






            new_cols = []
            if ab_set in ['signaling1', 'signaling2']:
                for channel in signaling_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_signaling[ab]+'_Nuc') #This file needs to be moved and imported
                for channel in signaling_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_signaling[ab]+'_Cyto') #This file needs to be moved and imported
            else:
                for channel in proteome_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_proteome[ab]+'_Nuc') #This file needs to be moved and imported
                for channel in proteome_channels: #Change this to detect correct channel lsit
                    ab = (channel,'row'+row)
                    new_cols.append(ab_dict_proteome[ab]+'_Cyto') #This file needs to be moved and imported

            new_cols = df.columns + optDF.columns
            df.columns = new_cols




ab_dict_proteome={('DAPI-0001','rowB'):'DNA-1',
('DAPI-0002','rowB'):'DNA-2',
('DAPI-0003','rowB'):'DNA-3',
('DAPI-0004','rowB'):'DNA-4',
('FITC-0001','rowB'):'MLC2v',
('FITC-0002','rowB'):'PCNA',
('FITC-0003','rowB'):'Ki67',
('FITC-0004','rowB'):'p-STAT5',
('Cy3-0001','rowB'):'MCM6',
('Cy3-0002','rowB'):'Oct4a',
('Cy3-0003','rowB'):'p-Rb',
('Cy3-0004','rowB'):'p-Aurora',
('Cy5-0001','rowB'):'a-ACTININ',
('Cy5-0002','rowB'):'Mlc2a',
('Cy5-0003','rowB'):'Cox4',
('Cy5-0004','rowB'):'mTor',
('DAPI-0001','rowC'):'DNA-1',
('DAPI-0002','rowC'):'DNA-2',
('DAPI-0003','rowC'):'DNA-3',
('DAPI-0004','rowC'):'DNA-4',
('FITC-0001','rowC'):'MLC2v',
('FITC-0002','rowC'):'PCNA',
('FITC-0003','rowC'):'Ki67',
('FITC-0004','rowC'):'p-STAT5',
('Cy3-0001','rowC'):'MCM6',
('Cy3-0002','rowC'):'Oct4a',
('Cy3-0003','rowC'):'p-Rb',
('Cy3-0004','rowC'):'p-Aurora',
('Cy5-0001','rowC'):'a-ACTININ',
('Cy5-0002','rowC'):'Mlc2a',
('Cy5-0003','rowC'):'Cox4',
('Cy5-0004','rowC'):'mTor',
('DAPI-0001','rowD'):'DNA-1',
('DAPI-0002','rowD'):'DNA-2',
('DAPI-0003','rowD'):'DNA-3',
('DAPI-0004','rowD'):'DNA-4',
('FITC-0001','rowD'):'MLC2v',
('FITC-0002','rowD'):'EGFR',
('FITC-0003','rowD'):'STAT3',
('FITC-0004','rowD'):'p-PKM2',
('Cy3-0001','rowD'):'HCN4',
('Cy3-0002','rowD'):'LC3',
('Cy3-0003','rowD'):'B-actin',
('Cy3-0004','rowD'):'PDGFRa',
('Cy5-0001','rowD'):'Flt1',
('Cy5-0002','rowD'):'Mlc2a',
('Cy5-0003','rowD'):'B-catenin',
('Cy5-0004','rowD'):'HIF1a',
('DAPI-0001','rowE'):'DNA-1',
('DAPI-0002','rowE'):'DNA-2',
('DAPI-0003','rowE'):'DNA-3',
('DAPI-0004','rowE'):'DNA-4',
('FITC-0001','rowE'):'MLC2v',
('FITC-0002','rowE'):'EGFR',
('FITC-0003','rowE'):'STAT3',
('FITC-0004','rowE'):'p-PKM2',
('Cy3-0001','rowE'):'HCN4',
('Cy3-0002','rowE'):'LC3',
('Cy3-0003','rowE'):'B-actin',
('Cy3-0004','rowE'):'PDGFRa',
('Cy5-0001','rowE'):'Flt1',
('Cy5-0002','rowE'):'Mlc2a',
('Cy5-0003','rowE'):'B-catenin',
('Cy5-0004','rowE'):'HIF1a',
('DAPI-0001','rowF'):'DNA-1',
('DAPI-0002','rowF'):'DNA-2',
('DAPI-0003','rowF'):'DNA-3',
('DAPI-0004','rowF'):'DNA-4',
('FITC-0001','rowF'):'MLC2v',
('FITC-0002','rowF'):'p53',
('FITC-0003','rowF'):'c-Jun',
('FITC-0004','rowF'):'Bax',
('Cy3-0001','rowF'):'p-MEK',
('Cy3-0002','rowF'):'pRb',
('Cy3-0003','rowF'):'c-Myc',
('Cy3-0004','rowF'):'B-tubulin',
('Cy5-0001','rowF'):'CNX43',
('Cy5-0002','rowF'):'p21',
('Cy5-0003','rowF'):'p27',
('Cy5-0004','rowF'):'Bcl2',
('DAPI-0001','rowG'):'DNA-1',
('DAPI-0002','rowG'):'DNA-2',
('DAPI-0003','rowG'):'DNA-3',
('DAPI-0004','rowG'):'DNA-4',
('FITC-0001','rowG'):'MLC2v',
('FITC-0002','rowG'):'p53',
('FITC-0003','rowG'):'c-Jun',
('FITC-0004','rowG'):'Bax',
('Cy3-0001','rowG'):'p-MEK',
('Cy3-0002','rowG'):'pRb',
('Cy3-0003','rowG'):'c-Myc',
('Cy3-0004','rowG'):'B-tubulin',
('Cy5-0001','rowG'):'CNX43',
('Cy5-0002','rowG'):'p21',
('Cy5-0003','rowG'):'p27',
('Cy5-0004','rowG'):'Bcl2'}

