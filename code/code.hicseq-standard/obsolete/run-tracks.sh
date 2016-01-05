#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-tracks.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running tracks ============="
resources=2,40G
method=by-group
scripts-master-loop.sh $resources $method ./code/hicseq-tracks.tcsh results/tracks "params/params.*.tcsh" "filter/results"


