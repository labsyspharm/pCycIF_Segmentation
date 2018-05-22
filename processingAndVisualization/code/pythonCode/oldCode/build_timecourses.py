#Folder for each drug
#One file per dose
#Time courses of all observables
#Can keep observables horizontal for now, would prefer vertical
#Nuc and Cyt in file but separate
#Abs in same file or separate?



import pandas as pd
import argparse
import glob
import sys 
from ab_dict_signaling_file import ab_dict_signaling
from ab_dict_proteome_file import ab_dict_proteome
from drug_code_dict import drug_code

#Allowed inputs
drugchoices = ['Afatinib','Crizotinib','Dasatinib','DMSO','Gefitinib','Lapatinib','Nilotinib','Ponatinib','Sorafenib','Sunitinib']
dosechoices = ['1','3.16','10']
timechoices = ['1','2','4','24','48','72','120']
abchoices = ['signaling1','signaling2','proteome1','proteome2','proteome3']

parser = argparse.ArgumentParser(description='Write')
parser.add_argument('--drugs', type=str, nargs='+', help='<Required> List of space-separated drug names. Allowed values are '+', '.join(drugchoices), metavar='drug names', choices=drugchoices, required=True)    
parser.add_argument('--dose', nargs='+', help='<Required> List of doses. Allowed values are '+', '.join(dosechoices), metavar='doses', choices=dosechoices, required=False)
parser.add_argument('--time', nargs='+', help='<Required> List of time points. Allowed values are '+', '.join(timechoices), metavar='times', choices=timechoices, required=False)
parser.add_argument('--abs', nargs='+', help='<Required> List of antibody sets. Allowed values are '+', '.join(abchoices), metavar='antibodies', choices=abchoices, required=False)
args = parser.parse_args()

drug_init = args.drugs
dose_init = args.dose
time_pts = args.time

dose_dict = {('10','signaling1'):'B',('10','signaling2'):'E',('3.16','signaling1'):'C',('3.16','signaling2'):'F',('3.16','proteome1'):'B',('3.16','proteome2'):'D',('3.16','proteome3'):'F',('1','signaling1'):'D',('1','signaling2'):'H',('1','proteome1'):'C',('1','proteome2'):'E',('1','proteome3'):'G'}
time_dict = {'1':'plate1/','2':'plate2/','4':'plate3/','24':'plate4/','48':'plate6/','72':'plate7/','96':'plate8/'}


plates_init = []
for t in time_pts:
    plates_init.append(time_dict[t])
ab_sets =  args.abs

drugs=[]
for d in drug_init:
    drugs.append(d+'_data/')


rows_sig_init = []
rows_prot_init = []
for dose in dose_init:
    for ab in ab_sets:
        if 'sig' in ab:
            rows_sig_init.append(dose_dict[(dose,ab)])
        else:
            if dose != '10':
                rows_prot_init.append(dose_dict[(dose,ab)])

#Function to filter rows, plates, based on ab_set
#Call withing main function
def filter_rows(ab_set,rows_sig,rows_prot):
    if 'sig' in ab_set:
        rows = rows_sig
    else:
        rows = rows_prot

    if ab_set == 'signaling1':
        rowlist = ['B','C','D']
    elif ab_set == 'signaling2':
        rowlist = ['E','F','G']
    elif ab_set == 'proteome1':
        rowlist = ['B','C']
    elif ab_set == 'proteome2':
        rowlist = ['D','E']
    elif ab_set == 'proteome3':
        rowlist = ['F','G']

    rows_to_use = list(set(rowlist) & set(rows))
    return rows_to_use

def filter_plates(ab_set,plates):
    if ab_set == 'signaling1' or ab_set == 'signaling2':
        platelist = ['plate1/','plate2/','plate3/','plate4/']
    elif ab_set == 'proteome1' or ab_set == 'proteome2' or ab_set == 'proteome3':
        platelist = ['plate6/','plate7/','plate8/']

    plates_to_use = list(set(platelist) & set(plates))
    return plates_to_use


def make_file_one_abset(ab_set):
#    print('input rows are %s' % rows_init)
    plates = filter_plates(ab_set,plates_init)
    rows = filter_rows(ab_set,rows_sig_init,rows_prot_init)
#    print('filtered rows are %s' % rows)
    if rows != [] and plates != []:


        #loop through plates
        for plate in plates:
            #loop through both subcellular locations    
            for loc in ['Nuc','Cyto']:
                file = 1
                #find all data files for a given plate and location, within folder for input drug

#Drug only accounted for in file name, need to loop through
                for drug in drugs:

                    for fn in sorted(glob.glob(drug+plate+'*'+loc+'*.txt')):
                        #limit to files that much input doses
           
                        for idx, row in enumerate(rows):
#                            print(rows)
#                            print(idx)
#                            print(ab_set)
                            if row == fn.split('-')[2][0]:  #CHanged row to list, to account for multiple input doses
        #                        print(fn)

                                df = pd.read_table(fn)
                                #write new dataframe with only row label and mean intensity, this is the only data we need for now
                                newdf = df.filter(['Label','Mean'])
                                #create list of all row labels
                                labellist = list(newdf['Label'])
                                secondarylist = []
                                cellnolist = []
                                #From each row label extract the cell number, and the secondary ab with cycle number and row letter
                                for label in labellist:
                                    secondarylist.append((label.split(':')[2],'row'+label[0])) ########HERE#################
                                    cellnolist.append(label.split(':')[1].split('-')[0])
                                observables = []
                                for ab in secondarylist:
                                    observables.append(ab_dict_signaling[ab])  #######Careful here
                                #Add columns in dataframe with cell numbers and biological observable
                                newdf['Cell_No'] = cellnolist
                                newdf['Observable'] = observables

                                #Want final dataframe to have one column for each biological observable, with values being mean intensity
                                newcols = list(set(observables))

                                #Build dataframe for each field file
                                i = 1
                                finaldf = pd.DataFrame()
                                for ob in newcols:
                                    #For first observable need to also record cell number
                                    if i == 1:
                                        finaldf['Cell_No'] = list(newdf.loc[newdf['Observable']==ob,'Cell_No'])
                                        finaldf[ob] = list(newdf.loc[newdf['Observable']==ob,'Mean'])
                                        i = i+1
                                    else:
                                        finaldf[ob] = list(newdf.loc[newdf['Observable']==ob,'Mean'])

                                finaldf['Plate'] = plate[5] 
                                finaldf['Dose'] = dose_init[idx] 
                                col = fn.split('-')[2][1:3]
                                finaldf['Column'] = col
                                finaldf['Field'] = fn[-5]
                                code = drug_code[(col,dose_init[idx])]   
                                finaldf['Dose_Code'] = code
                                #Stitch all fields into one dataframe
                                if file == 1:
                                    file = 2
                                    finaldf2 = finaldf
                                else:
                                    finaldf2 = pd.concat([finaldf2,finaldf],axis=0,ignore_index=True)

                #Combining nuc and cyt
                new_cols = []
                cols = list(finaldf2.columns)
                for col in cols:
                    new_cols.append(col+'_'+loc)

                if loc == 'Nuc':
                    combined1 = finaldf2
                    combined1.columns = new_cols #Check order doesn't change
                else:
                    combined_df = pd.DataFrame()
                    combined2 = finaldf2
                    combined2.columns = new_cols #Check order doesn't change
                    combined_df = pd.concat([combined1,combined2],axis=1,ignore_index=True)
                    combined_df.columns = list(combined1.columns) + list(combined2.columns)
                    combined_df = combined_df[sorted(combined_df.columns)]


            #at end of each plate 1-4
            #Print to final location
            print(plate)
            if plate == plates[0]: #Problem here. Can't just be one or six, needs to be first in list
                allplates_df = combined_df
    #        elif allplates_df: #is not None:
            else:
                try: 
                    allplates_df
                except UnboundLocalError:
                    pass
                else:
                    allplates_df = pd.concat([allplates_df,combined_df])
        ##filename = sys.argv[1]+'ab_sets_dosedependent/signaling1_%s.csv' % sys.argv[2]  #here #currently drug/intermediatefolder/abset_dose.csv
        #filename = 'testing/'+drug_init+'_dose'+dose+'_time'+time_pts[0]+ab_set+'_.csv'
        #allplates_df.to_csv(filename)

        drugstoprint = '_'.join(drug_init)
        timestoprint = '_'.join(time_pts)
        dosestoprint = '_'.join(dose_init)
        abstoprint = ab_set
#        filename = 'output/'+'drugs_%s_timepts_%s_doses_%s_absets_%s.csv' % (drugstoprint,timestoprint,dosestoprint,abstoprint)
        filename = '../testing_output/'+'drugs_%s_timepts_%s_doses_%s_absets_%s.csv' % (drugstoprint,timestoprint,dosestoprint,abstoprint)
    #    if all_plates_df is not None:
    #        allplates_df.to_csv(filename)
        try: 
            allplates_df
        except UnboundLocalError:
            pass
        else:
            allplates_df.to_csv(filename)

#######################################################################################3


for ab_set in ab_sets:
    make_file_one_abset(ab_set)

