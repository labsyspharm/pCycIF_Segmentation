#!/bin/bash
#BSUB -n 1                #Number of corse
#BSUB -W 12:00            #Wall time
#BSUB -J resegmentImages         #Job name
#BSUB -N -o %resegmentImages.out     #lsf output file
#BSUB -N -e %resegmentImages.err    #lsf error file
#BSUB -q short            #queue
#BSUB -w "done(copyData)"

source ~/venv3/bin/activate
for folder in ProcessedDataP1 ProcessedDataP2 ProcessedDataP3 ProcessedDataP4 ProcessedDataP5 ProcessedDataP6 ProcessedDataP7 ProcessedDataP8;
do rm /n/scratch2/rps21/ResegmentedPlates/$folder/*mask* ;
rm /n/scratch2/rps21/ResegmentedPlates/$folder/*.txt;
python resegmentImages_dapi.py --plate $folder;
python normalizeResegmentedImages_dapi.py --plate $folder
python logResegmentedImages_dapi.py --plate $folder
python logNormalizeResegmentedImages_dapi.py --plate $folder;
done;
