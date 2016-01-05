#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-filter-stats.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running filter-stats ============="
threads=1
method=by-branch
scripts-master-loop.sh $threads $method ./code/hicseq-filter-stats.tcsh results/filter-stats "params/params.*.tcsh" "filter/results"


