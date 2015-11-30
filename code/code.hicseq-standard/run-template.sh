#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-%%op%%.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running %%op%% ============="
threads=%%threads%%
method=%%method%%
scripts-master-loop.sh $threads $method ./code/hicseq-%%op%%.tcsh results/%%op%% "params/params.*.tcsh" "%%inpresults%%"


