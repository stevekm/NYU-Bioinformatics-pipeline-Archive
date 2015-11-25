#!/bin/bash
source ./code/code.main/custom-bashrc      # customize shell environment

##
## USAGE: run-align.sh
##

# process command-line inputs
if [ $# != 0 ]
then
  grep '^##' $0
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Aligning reads ============="
threads=16
scripts-master-loop.sh $threads by-sample ./code/hicseq-align.tcsh results/align "params/params.*.tcsh" inputs/fastq


