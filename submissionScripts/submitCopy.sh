#!/bin/bash
#BSUB -n 1                #Number of corse
#BSUB -W 2:00            #Wall time
#BSUB -J copyData         #Job name
#BSUB -N -o %copyData.out     #lsf output file
#BSUB -N -e %copyData.err    #lsf error file
#BSUB -q short            #queue
#BSUB -R "select[transfer]" 


for folder in ProcessedDataP1 ProcessedDataP2 ProcessedDataP3 ProcessedDataP4 ProcessedDataP5 ProcessedDataP6 ProcessedDataP7 ProcessedDataP8;
do cp -r /files/ImStor/sorger/data/IN\ Cell\ Analyzer\ 6000/Connor/June\ FDA\ CycIF\ Processed/Processed\ Plates/$folder /n/scratch2/rps21/ResegmentedPlates/$folder;
done;
