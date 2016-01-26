#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-matrix-filtered.tcsh [--dry-run]
##

#% This step generates Hi-C contact matrices, per chromosome and/or genome-wide according to specified parameters and for all desired resolutions. Read-pairs that passed the filtering step were used to create a contact matrix. Each element of the matrix corresponds to a pair of non-overlapping genomic bins. The value of each element corresponds to the number of paired reads overlapping the corresponding genomic binds and the size of each one of the bins corresponds to the resolution. 
#% TABLES: 
#% FIGURES: 

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit
endif

set opt = "$1"

# setup
set op = matrix-filtered
set inpdirs = "inpdirs/*"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 2,20G
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/hicseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "sample" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



