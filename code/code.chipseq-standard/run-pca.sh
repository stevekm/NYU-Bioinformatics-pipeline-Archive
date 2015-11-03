#!/bin/bash

##
## USAGE: run-pca.sh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# perform PCA
scripts-send2err "=== Performing PCA ============="
scripts-create-path results/
threads=1
op=pca
scripts-master-loop.sh $threads by-branch ./code/chipseq-$op.tcsh results/$op "params/params.*.tcsh" matrices/results


