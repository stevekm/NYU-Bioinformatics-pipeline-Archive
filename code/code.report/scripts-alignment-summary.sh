#!/bin/bash

##
## USAGE: scripts-alignment-summary.sh OutputDir ProjectDir
##

#% DESCRIPTION: This wrapper script will aggregate summary alignment statistics and create barplots

# make sure that the correct number of script arguments were used
if [ $# != 2 ]
then
  grep '^##' $0 # print out lines from the script that start with ##
  exit
fi


module unload r/3.0.2
module load r/3.2.0 # need this version of R, always!

# get the path to the R companion script
BarplotScript=$(readlink -m $(dirname $0)) # get the current script path
CurrentScript=$(basename $0)
BarplotScript=$BarplotScript/${CurrentScript/.sh/.R}

# get the script arguments
OutputDir=$1
ProjDir=$2

# make sure the OutputDir exists
mkdir -p $OutputDir

# get the path to the alignment pipeline summary object
AlignDB=$ProjDir/pipeline/align/results/.db/db.RData

# get the path to the alignment directory
AlignDir=$ProjDir/pipeline/align

# make sure to call Rscript directly here to get the correct version and args entries in the script!
Rscript $BarplotScript "$OutputDir" "$AlignDB" "$AlignDir"

