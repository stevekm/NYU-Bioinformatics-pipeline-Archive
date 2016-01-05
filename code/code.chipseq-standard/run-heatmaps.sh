#!/bin/bash
source ./code/code.main/custom-bashrc    # shell settings

##
## USAGE: run-heatmaps.sh
##


# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# perform PCA
scripts-send2err "=== Generating heatmaps ============="
scripts-create-path results/
threads=1
op=heatmaps
scripts-master-loop.sh $threads by-branch ./code/chipseq-$op.tcsh results/$op "params/params.*.tcsh" "inpdirs/matrices/results/matrices.*.nbins_25"


