#!/bin/bash
source ./code/code.main/custom-bashrc     # setup bash environment

##
## USAGE: run-pca.sh
##

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# perform PCA
scripts-send2err "=== Performing QC ============="
scripts-create-path results/
threads=4
scripts-master-loop.sh $threads by-group ./code/chipseq-fingerprint.tcsh results/fingerprint "params/params.*.tcsh" "inpdirs/*/results"


