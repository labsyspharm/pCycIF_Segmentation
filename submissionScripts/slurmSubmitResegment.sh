#!/bin/bash
#SBATCH --mail-type=ALL 	  #Sends email for options BEGIN, END, FAIL, ALL
#SBATCH -c 1			  #Number of cores
#SBATCH -N 1                	  #Number of nodes
#SBATCH -t 0-2:00            	  #Runtime in D-HH:MM format
#SBATCH -J resegmentImages        #Job name
#SBATCH -o %resegmentImages.out	  #stdout file
#SBATCH -e %resegmentImages.err   #stderr file
#SBATCH -p priority            	  #Partition


/home/rps21/pCycIF_Segmentation/run_runSegmentationO2.sh /home/rps21/MCR/v901/ externalParams.txt





##source ~/venv36/bin/activate
##for folder in ProcessedDataP2 ProcessedDataP3 ProcessedDataP4 ProcessedDataP6 ProcessedDataP7 ProcessedDataP8;
###for folder in ProcessedDataP1
##do rm /n/scratch2/rps21/ResegmentedPlates/$folder/_* ;
##python ../resegmentImages.py --plate $folder;
###python ../processResegmentedImages.py --plate $folder;
##done;
