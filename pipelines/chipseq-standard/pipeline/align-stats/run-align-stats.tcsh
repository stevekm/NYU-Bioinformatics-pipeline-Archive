#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-align-stats.tcsh [--dry-run]
##

#TITLE: Alignment Summary Stats
#DESCRIPTION: Create barplots to show the number of reads that are aligned versus unaligned. 
#FIGURE: alignment_barplots.pdf
#SHEET: alignment_stats_extended.csv

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set opt = "$1"

# setup
set op = align-stats
set inpdirs = "inpdirs/*"
set filter = ""
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 1
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v -F "$filter" "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "*" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



