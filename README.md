# Dockerized app for segmentation of plate-based Cyclic Immunofluorescence (pCycIF) imaging data

## Segmentation Methodology
Add description

## Installation
Download the docker image using the command `docker pull rps21/pcycif_segmentation`
Build the docker image by `docker build -t rps21/pcycif_segmentation .`

## Usage 
Create an input folder which contains one folder for each plate that was imaged. Each of these folders must contain the images to be segmented
Example:
/localPath/Plate1/image1.tif
/localPath/Plate1/image2.tif
/localPath/Plate2/image1.tif

Create a configuration folder, containing the file pcycif_segmentation.yml
All configuration parameters must be set by the user and are explained below. Example configuration file can be found at 
https://github.com/sorgerlab/pCycIF_Segmentation/blob/master/dockerBuild/config/cycif_segmentation.yml

Run the docker container, indicating paths to the input images, configuration file, and desired output location
`docker run -v /local/path/to/config/:/config -v /local/path/to/input/:/input -v /local/path/to/output/:/output rps21/pcycif_segmentation`


# Application configuration and example 
A full example, as well as the segmentation source code, is contained in the associated github repository which can be cloned locally using:
git clone https://github.com/sorgerlab/pCycIF_Segmentation.git

## Input Data
Input data must be .tif image stacks, organized in folders by plate imaged. 
Example data can be found at: /localPath/pCycIF_Segmentation/dockerBuild/input

## Segmentation Output
Output folder must be specified, and can correspond to any already created local folder. Output files consist of...
Example output can be found at:/localPath/pCycIF_Segmentation/dockerBuild/output


## Configuration parameters
The following parameters must be set by the user based on their experimental design. 
Example configuration file can be found at: /localPath/pCycIF_Segmentation/dockerBuild/config

parallel: Boolean value (0/1) determining if segmentation will be run in parallel
Example: 0
NucMaskChan: Matrix with numbered index of first and last image of tiff stack that correspond to nuclear stains
Example: [1 9]
CytoMaskChan: Matrix with numbered index of first and last image of tiff stack that correspond to cytoplasmic stains
Example: [10 36]
numCycles: Number value corresponding to the number of cycles
Example: 9
row: String containing the lettered index corresponding to the first and last row of a 96-well plate to segment 
Example: 'CE'
col: Matrix with numbered index of first and last column  of a 96-well plate to segment
Example: [4 10]
saveFig: Boolean value (0/1) determining if you want to save Matlab .fig files 
Example: 1
cytoMethod: Cytoplasm segmentation method to use. 
Example: 'RF'
MedianIntensity: 1
saveMasks: Boolean value (0/1) determining if you want to save mask images
Example: 1
applyFFC: String defining flat feel correction to use. Options:'none', 'ffonly'
Example: 'ffonly'
useRFNuc: Boolean value (0/1) determining if you want to use... ?
Example: 1
segmentCytoplasm: 'segmentCytoplasm'
bleached: 0


