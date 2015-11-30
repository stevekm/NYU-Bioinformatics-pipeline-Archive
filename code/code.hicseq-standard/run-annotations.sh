#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-annotations.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running annotations ============="
threads=1
method=by-object
scripts-master-loop.sh $threads $method ./code/hicseq-annotations.tcsh results/annotations "params/params.*.tcsh" "matrix.filtered/results/matrix_filtered.res_40kb/."


