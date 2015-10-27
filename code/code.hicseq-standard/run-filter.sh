#!/bin/bash

##
## USAGE: run-filter.sh
##

if [ $# != 0 ]
then
  grep '^##' $0
  exit
fi

# set path
PATH=./code/code:$PATH

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Filtering alignments ============="
threads=2
scripts-master-loop.sh $threads ./code/hicseq-filter.tcsh results/filter "params/params.*.tcsh" align/results


