#!/bin/bash

##
## USAGE: run-peaks.sh
##

if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# set path
PATH=./code/code:$PATH

# generate peaks
scripts-send2err "=== Generating peaks ============="
scripts-create-path results/
sheet=inputs/sample-sheet.tsv
threads=1
jid=()
for params in $(ls -1 params/params.*); do
  # create output directory for each parameter file
  out=results/peaks.$(echo $params | sed 's/.*\///' | sed 's/^params\.//')
  scripts-send2err "- out = $out"
  scripts-create-path $out/

  # call peaks for each alignment directory
  for aln in $(ls -1d alignments/results/align.*); do
    sample=$(echo $aln | sed 's/.*\/align\.//')                                        # "treatment" sample name
    bam_sample=$aln/alignments.bam
    control=$(cat $sheet | grep -v '^#' | cut -f2,3 | grep "^$sample	" | cut -f2)     # control (Input/IgG) sample name
    bam_control= 
    if [ "$control" != 'n/a' ]; then
      bam_control=alignments/results/align.$control/alignments.bam
    fi
    outdir=$out/peaks.$sample
    jid+=( $(scripts-qsub-wrapper $threads ./code/peaks-call $outdir $params $bam_sample $bam_control) )
  done
done

# done
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



