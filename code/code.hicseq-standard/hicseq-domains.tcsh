#!/bin/tcsh

source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH SAMPLE
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set sample = $4

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# Run domains
set input_sample = ($branch/$sample)

if ($tool == armatus) then
  ./code/hicseq-domains-$tool.tcsh $outdir "$input_sample" $gamma
else
  scripts-send2err "Error: unknown domain caller tool $tool."
  exit 1
endif

