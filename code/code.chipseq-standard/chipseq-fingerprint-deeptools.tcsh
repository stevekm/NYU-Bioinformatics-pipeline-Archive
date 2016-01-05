#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: chipseq-deeptools-fingerprint.tcsh OUTDIR PARAM-FILE INPUT-BAM-FILES LABELS
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

# Input
set outdir = $1
set params = $2
set input_bams = ($3)
set labels = ($4)

# setup
scripts-create-path $outdir/
source $params

# determine number of slots
if (! $?NSLOTS) then
  set threads = 8
else
  set threads = $NSLOTS
endif

# Test chip quality
bamFingerprint --verbose -p $threads --bamfiles $input_bams --labels $labels --plotFile $outdir/chip-fingerprint.pdf



