#!/bin/bash

##
## USAGE: run-matrix.sh
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

# setup
sheet=inputs/sample-sheet.tsv
release=inputs/release
cat $release/../DNA/genome.bed | cols -t 0 1 2 0 > results/genome.bed

# matrix
scripts-send2err "=== Generating contact matrices ============="
threads=2
scripts-master-loop.sh $threads by-group ./code/hicseq-matrix.tcsh results/matrix "params/params.*.tcsh" filter/results


