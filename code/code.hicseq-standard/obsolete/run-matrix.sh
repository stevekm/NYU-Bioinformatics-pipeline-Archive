#!/bin/bash

##
## USAGE: run-matrix.sh
##

# shell settings
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]
then
  grep '^##' $0
  exit
fi

# create results directory
scripts-create-path results/

# matrix
scripts-send2err "=== Generating contact matrices ============="
threads=2,20G
scripts-master-loop.sh $threads by-sample ./code/hicseq-matrix.tcsh results/matrix-filtered "params/params.*.tcsh" filter/results


