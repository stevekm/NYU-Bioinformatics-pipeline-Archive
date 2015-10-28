#!/bin/bash

##
## USAGE: run-filter.sh
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

# filter
scripts-send2err "=== Filtering alignments ============="
threads=2
scripts-master-loop.sh $threads ./code/hicseq-filter.tcsh results/filter "params/params.*.tcsh" align/results


