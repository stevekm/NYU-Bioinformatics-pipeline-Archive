#!/bin/bash
source ./code/code.main/custom-bashrc        # shell settings (must be included in all scripts)

##
## USAGE: run-matrices.sh 
##

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# generate matrices
scripts-send2err "=== Generating matrices ============="
scripts-create-path results/
sheet=inputs/sample-sheet.tsv
threads=1
op=matrices
scripts-master-loop.sh $threads by-sample ./code/chipseq-$op.tcsh results/$op "params/params.*.tcsh" "inpdirs/peaks/results/by-sample"







