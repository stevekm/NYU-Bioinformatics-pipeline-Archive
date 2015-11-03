#!/bin/bash

##
## USAGE: run-peaktable.sh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# generate matrices
scripts-send2err "=== Generating peak tables ============="
scripts-create-path results/
threads=1
op=peaktable
sh=sh
scripts-master-loop.sh $threads by-branch ./code/chipseq-$op.$sh results/$op "params/params.*.$sh" peaks/results




