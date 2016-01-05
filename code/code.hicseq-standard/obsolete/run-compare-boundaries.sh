#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-compare-boundaries.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running compare-boundaries ============="
threads=1
method=by-branch
scripts-master-loop.sh $threads $method ./code/hicseq-compare-boundaries.tcsh results/compare-boundaries "params/params.*.tcsh" "domains/results"


