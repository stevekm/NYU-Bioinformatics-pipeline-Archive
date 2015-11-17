#!/bin/bash
source ./code/code.main/custom-bashrc     # setup bash environment

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
scripts-send2err "=== Performing QC ============="
scripts-create-path results/
threads=4
scripts-master-loop.sh $threads by-group ./code/chipseq-fingerprint.sh results/fingerprint "params/params.*.sh" alignments/results
#scripts-master-loop.sh $threads by-branch ./code/chipseq-pca.tcsh results/pca "params/params.*.tcsh" matrices/results


