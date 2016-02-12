#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-boundary-scores-pca.tcsh [--dry-run]
##

#% This step generates PCA plots on boundary scores. 
#% TABLES: 
#% FIGURES: *.pdf

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit
endif

set opt = "$1"

# setup
set op = boundary-scores-pca
set inpdirs = "inpdirs/*"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 1
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/hicseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "genome" "*" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



