#!/bin/bash
source ./code/code.main/custom-bashrc             # shell settings (must be included in all scripts)

##
## USAGE: run-peaks.sh
##

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0 | scripts-send2err
  exit 1
fi

# create results directory
scripts-create-path results/

# generate peaks
scripts-send2err "=== Generating peaks ============="
threads=1
op=peaks
for method in "by-sample" "by-group"; do
  scripts-master-loop.sh $threads $method ./code/chipseq-$op.tcsh results/$method/$op "params/params.*.tcsh" "inpdirs/*/results"
done




