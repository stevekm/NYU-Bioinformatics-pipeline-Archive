#!/bin/bash
source ./code/code.main/custom-bashrc      # shell settings (must be included in all scripts)

##
## USAGE: run-peaktable.sh
##

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# generate matrices
scripts-send2err "=== Generating peak tables ============="
scripts-create-path results/
threads=1
op=peaktable
scripts-master-loop.sh $threads by-branch ./code/chipseq-$op.tcsh results/$op "params/params.*.tcsh" "inpdirs/*/results"




