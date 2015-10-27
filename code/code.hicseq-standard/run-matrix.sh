#!/bin/bash

##
## USAGE: run-matrix.sh
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

# setup
sheet=inputs/sample-sheet.tsv
groups=$(cat $sheet | grep -v '^#' | cut -f4 | sort -u)
release=inputs/release
cat $release/../DNA/genome.bed | cols -t 0 1 2 0 > results/genome.bed

# matrix
scripts-send2err "=== Generating contact matrices ============="
threads=2
scripts-master-loop-grouped.sh $threads ./code/hicseq-matrix.tcsh results/matrix "params/params.*.tcsh" filter/results "$groups"


