#!/bin/bash

##
## USAGE: run-peaks.sh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# create results directory
scripts-create-path results/

# generate peaks
scripts-send2err "=== Generating peaks ============="
threads=1
op=peaks
for method in "by-sample" "by-group"; do
  scripts-master-loop.sh $threads $method ./code/chipseq-$op.tcsh results/$method/$op "params/params.*.tcsh" alignments/results
done




