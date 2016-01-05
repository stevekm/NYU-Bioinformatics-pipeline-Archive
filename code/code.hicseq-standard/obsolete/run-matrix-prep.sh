#!/bin/bash

##
## USAGE: run-matrix-prep.sh
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
scripts-send2err "=== Preprocess matrices ============="
threads=1,40G      # TODO: automatically determine memory usage
op=matrix-prep
scripts-master-loop.sh $threads by-object ./code/hicseq-$op.tcsh results/$op "params/params.*.tcsh" "matrix-filtered/results/*_40kb"


