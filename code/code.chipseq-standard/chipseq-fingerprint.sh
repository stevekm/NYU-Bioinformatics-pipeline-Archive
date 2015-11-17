#!/bin/bash
source ./code/code.main/custom-bashrc     # setup bash environment

##
## USAGE: chipseq-fingerprint.sh OUTPUT-DIR PARAMETER-SCRIPT ALIGNMENT-BRANCH SAMPLES
##

# process command-line inputs
if (($# != 4)); then
  grep '^##' $0
  exit
fi

# inputs
outdir=$1
params=$2
aln_branch=$3
samples=$4

# create output path
scripts-create-path $outdir/

# determine input files
treatment_bam=$(echo $samples | tr ' ' '\n' | awk -v d=$aln_branch '{print d"/"$0"/alignments.bam"}')         # TODO: create function in query-sample-sheet.tcsh
control_samples=$(./code/query-sample-sheet.tcsh control "$samples" | sort -u) 
if [ "$control_samples" == "" ]; then
  control_bam=
else
  control_bam=$(echo $control_samples | tr ' ' '\n' | awk -v d=$aln_branch '{print d"/"$0"/alignments.bam"}')   # TODO: create function in query-sample-sheet.tcsh
fi
bam=$(echo $treatment_bam $control_bam)
labels=$(echo $samples $control_samples)

# run deeptools wrapper
./code/chipseq-fingerprint-deeptools.sh $outdir $params "$bam" "$labels"



