#!/bin/bash

##
## USAGE: run-ic.sh
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
resources=1,80G
inpdirs=$(ls -1d matrix-filtered/results/* | grep -v '\.maxd_')                                                             # do not run IC on distance-restricted matrices
scripts-master-loop.sh $resources by-object ./code/hicseq-ic.tcsh results/matrix-ic "params/params.*.tcsh" "$inpdirs"


