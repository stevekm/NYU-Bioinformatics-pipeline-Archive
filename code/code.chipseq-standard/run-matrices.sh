#!/bin/bash

##
## USAGE: run-matrices.sh 
##

if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# set path
PATH=./code/code:$PATH

# generate matrices
scripts-send2err "=== Generating matrices ============="
scripts-create-path results/
threads=1
jid=()
# loop over all matrix parameter settings
for params in $(ls -1 params/params.*.sh); do
 
  # loop over all peak parameter settings
  for peak_params in $(find peaks/results/ -name 'peaks.bed' | sed 's/\/peaks\.[^/]\+\/peaks.bed//' | sort -u); do
    # create output directory
    out=results/matrices.$(echo $params | sed 's/.*\///' | sed 's/^params\.//' | sed 's/\.sh$//')/$(echo $peak_params | sed 's/peaks\/results\///')
    ref_regions=( $(ls -1 $peak_params/peaks.*/peaks.bed) )
    scripts-send2err "-- out = $out"

    # create matrices
    for aln in $(ls -1d alignments/results/align.*); do
      sample=$(echo $aln | sed 's/.*\/align\.//')
      jid+=( $(scripts-qsub-wrapper $threads ./code/create-matrix.sh $out/matrix.$sample $params $aln/alignments.bam "${ref_regions[*]}") )
    done
  done
done

# wait for all jobs to finish
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"

# done
scripts-send2err "Done."






