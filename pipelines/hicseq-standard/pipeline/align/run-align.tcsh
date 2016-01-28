#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-align.tcsh [--dry-run]
##

#% This step performs alignment of paired-end Hi-C reads. 
#% TABLES: 
#% FIGURES: 

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit
endif

set opt = "$1"

# setup
set op = align
set inpdirs = "inputs"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 16
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/hicseq-$op.tcsh"

# symlink fastq to results (required for first step of the pipeline)
if (! -e inputs/$results) then
  (cd inputs; ln -s fastq $results)
endif

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "sample" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



