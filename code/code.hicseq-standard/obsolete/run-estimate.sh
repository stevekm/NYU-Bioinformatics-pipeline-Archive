#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

##
## USAGE: run-estimate.sh
##

# shell settings
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# create results directory
scripts-create-path results/

# matrix
scripts-send2err "=== Estimating matrix using fused lasso ============="
resources=1
inpdirs=$(find matrix-*/results -name 'matrix-*.res_40kb')
scripts-master-loop.sh $resources by-object ./code/hicseq-estimate.tcsh results/matrix-estimated "params/params.*.tcsh" "$inpdirs"

