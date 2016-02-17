#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-align.tcsh [--dry-run]
##

# ~~~ Entries for auto-report ~~~ #
#TITLE: Alignment
#DESCRIPTION: This step takes fastq files as input for each sample, and perform alignment using different aligners according to the specified parameter files.
#FIGURE:
#SHEET:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

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
set resources = 8
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# symlink fastq to results (required for first step of the pipeline)
if (! -e inputs/$results) then
  (cd inputs; ln -s fastq $results)
endif

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "sample" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run
