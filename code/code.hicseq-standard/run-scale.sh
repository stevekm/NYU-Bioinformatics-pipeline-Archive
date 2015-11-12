#!/bin/bash

##
## USAGE: run-scale.sh
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
scripts-send2err "=== Generating scaled matrices ============="
threads=1
scripts-master-loop.sh $threads by-object ./code/hicseq-scale.tcsh results/matrix_scaled "params/params.*.tcsh" matrix.filtered/results


