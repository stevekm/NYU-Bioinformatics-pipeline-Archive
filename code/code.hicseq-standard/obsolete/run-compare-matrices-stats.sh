#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-compare-matrices-stats.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# run
op=compare-matrices-stats
scripts-send2err "=== Running $op ============="
threads=1
method=by-branch
scripts-master-loop.sh $threads $method ./code/hicseq-$op.tcsh results/$op "params/params.*.tcsh" "compare-matrices/results"


