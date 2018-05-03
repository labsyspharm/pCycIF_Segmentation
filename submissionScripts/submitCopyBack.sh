#!/bin/bash
#BSUB -n 1                #Number of corse
#BSUB -W 4:00            #Wall time
#BSUB -J copyBackFiles         #Job name
#BSUB -N -o %copyBackFiles.out     #lsf output file
#BSUB -N -e %copyBackFiles.err    #lsf error file
#BSUB -q short            #queue
#BSUB -w "done(removeFiles)"
#BSUB -R "select[transfer]"

for folder in ProcessedDataP1 ProcessedDataP2 ProcessedDataP3 ProcessedDataP4 ProcessedDataP5 ProcessedDataP6 ProcessedDataP7 ProcessedDataP8;
cp -r /n/scratch2/rps21/ResegmentedPlates/$folder/ /files/ImStor/sorger/data/IN\ Cell\ Analyzer\ 6000/Connor/June\ FDA\ CycIF\ Processed/ResegmentedPlates/$folder;
