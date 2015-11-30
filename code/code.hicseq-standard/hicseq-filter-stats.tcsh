#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-filter-stats.tcsh OUTPUT-DIR PARAM-SCRIPT FILTER-BRANCH
##

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# Plot barplots
set input_samples = ($branch/*)
Rscript ./code/hicseq-filter-stats.r $outdir "$input_samples"


