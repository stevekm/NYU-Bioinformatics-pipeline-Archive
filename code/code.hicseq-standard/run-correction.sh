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
scripts-send2err "=== Generating corrected matrix ============="
threads=1
scripts-master-loop.sh $threads by-object ./code/hicseq-ic.tcsh results/estimation "params/params.*.tcsh" matrix/results


