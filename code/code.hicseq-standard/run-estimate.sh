#!/bin/bash

##
## USAGE: run-estimate.sh
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
scripts-send2err "=== Estimating matrix using fused lasso ============="
resources=1
inpdirs=$(find matrix-*/ -name 'matrix-*.res_*kb')
scripts-master-loop.sh $resources by-object ./code/hicseq-estimate.tcsh results/matrix-estimated "params/params.prep_log2.fused1Dsymm.tcsh" "$inpdirs"
#scripts-master-loop.sh $resources by-object ./code/hicseq-estimate.tcsh results/matrix-estimated "params/params.*.tcsh" "matrix-*/results"

