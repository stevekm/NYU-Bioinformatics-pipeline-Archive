#!/bin/bash

##
## USAGE: run-peaks.sh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

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
    outdir=$out/peaks.$sample
    if [ $(cat $sheet | cut -f2 | grep -v '^#' | grep -c "^$sample$") == 0 ]; then
      scripts-send2err "Warning: sample $sample is not in the sample sheet, removing..."
      rm -rf $outdir
      continue
    fi
    bam_sample=$aln/alignments.bam
    control=$(cat $sheet | grep -v '^#' | cut -f2,3 | grep "^$sample	" | cut -f2)     # control (Input/IgG) sample name
    bam_control= 
    if [ "$control" != 'n/a' ]; then
      bam_control=alignments/results/align.$control/alignments.bam
    fi
    jid+=( $(scripts-qsub-wrapper $threads ./code/peaks-call $outdir $params $bam_sample $bam_control) )
  done
done

# done
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



