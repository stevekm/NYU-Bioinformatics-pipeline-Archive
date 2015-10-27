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

# align
scripts-send2err "=== Filtering alignments ============="
threads=2
jid=()
# loop over all parameter files
for filter_params in $(ls -1 params/params.*.tcsh); do
  pname=$(echo $filter_params | sed 's/.*\///' | cut -d'.' -f2)

  # loop over all alignments
  inpdir=align/results
  for aln in $(cd $inpdir; find . -name 'job.sh' | sed 's/\/job.sh$//' | sed 's/^.\///'); do
    outdir=results/filter.$pname/$aln
    jid+=( $(scripts-qsub-wrapper $threads ./code/hicseq-filter.tcsh $outdir $filter_params $inpdir/$aln/alignments.bam) )
  done

done
       
# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



