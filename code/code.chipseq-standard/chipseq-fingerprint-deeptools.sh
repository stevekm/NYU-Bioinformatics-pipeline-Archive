#!/bin/bash
source ./code/code.main/custom-bashrc     # setup bash environment

##
## USAGE: chipseq-deeptools-fingerprint.sh OUTDIR PARAM-FILE INPUT-BAM-FILES LABELS
##

# process command-line inputs
if (($# != 4)); then
  grep '^##' $0
  exit
fi

# Input
outdir=$1
params=$2
input_bams=$3
labels=$4

# setup
scripts-create-path $outdir/
source $params

if [ $?NSLOTS ]; then
  threads=$NSLOTS
else
  threads=1
fi

# Test chip quality
bamFingerprint --verbose -p $threads --bamfiles $input_bams --labels $labels --plotFile $outdir/chip-fingerprint.pdf



