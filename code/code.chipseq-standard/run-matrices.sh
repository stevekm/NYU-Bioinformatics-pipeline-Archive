#!/bin/bash

##
## USAGE: run-matrices.sh 
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

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
sh=sh
scripts-master-loop.sh $threads by-sample ./code/chipseq-$op.$sh results/$op "params/params.*.$sh" peaks/results/by-sample







