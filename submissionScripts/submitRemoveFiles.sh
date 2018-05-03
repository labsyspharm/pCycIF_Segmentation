#!/bin/bash
#BSUB -n 1                #Number of corse
#BSUB -W 1:00            #Wall time
#BSUB -J removeFiles         #Job name
#BSUB -N -o %removeFiles.out     #lsf output file
#BSUB -N -e %removeFiles.err    #lsf error file
#BSUB -q short            #queue
#BSUB -w "done(resegmentImages)"


shopt -s extglob
for folder in ProcessedDataP1 ProcessedDataP2 ProcessedDataP3 ProcessedDataP4 ProcessedDataP5 ProcessedDataP6 ProcessedDataP7 ProcessedDataP8;
rm /n/scratch2/rps21/ResegmentedPlates/$folder/!(_*) 




