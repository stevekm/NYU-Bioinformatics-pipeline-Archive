#!/bin/bash

##
## USAGE: scripts-master-loop-grouped.sh THREADS METHOD SCRIPT OUTDIR-PREFIX PARAM-SCRIPTS INPUT-DIR
##
##   THREADS         number of threads
##   METHOD          script will be run using one of these methods: by-object, by-sample, by-group, by-branch 
##   SCRIPT          script to be run for each combination of parameters and input branch, e.g. chipseq-peaks.tcsh
##   PARAM-SCRIPTS   array of parameter scripts
##   INPUT-DIR       directory of computation tree to be used as input
##

if (( ($# < 6) || ($# > 6) )); then
  grep '^##' $0
  exit
fi

threads=$1
method=$2
operation=$3
outpref=$4
paramset=($5)
inpdir=$6

# setup
sheet=inputs/sample-sheet.tsv
samples=( $(cat $sheet | grep -v '^#' | cut -f1) )
groups=( $(cat $sheet | grep -v '^#' | cut -f2 | tr '|' '\n' | sort -u) )

# determine input branches
branches=( $(cd $inpdir; find . -name 'job.sh' | sed 's/\/job.sh$//' | sed 's/\/[^/]\+$//' | sed 's/^\.\///' | sort -u) )
if [ ${#branches[*]} == 0 ]; then
  branches=( '.' )
fi

# start
jid=()
#
# loop over all parameter files
#
for p in ${paramset[*]}; do
  pname=$(echo $p | sed 's/.*\///' | sed 's/^params\.//' | sed 's/\.sh$//' | sed 's/\.tcsh$//')

  #
  # loop over all branches of input tree (one level before the leaves)
  #
  for branch in ${branches[*]}; do

    if [ $method == "by-sample" ]; then
      #
      # loop over all samples
      #
      for s in ${samples[*]}; do
        outdir=$outpref.$pname/$branch/$s
        jid+=( $(scripts-qsub-wrapper $threads $operation $outdir $p $inpdir/$branch $s) )
      done

    elif [ $method == "by-group" ]; then
      #
      # loop over all groups
      #
      for g in ${groups[*]}; do
        outdir=$outpref.$pname/$branch/$g
        gsamples=( $(cat $sheet | grep -v '^#' | tr '|' ' ' | key_expand | awk -v g=$g '$2==g' | cut -f1) )
        jid+=( $(scripts-qsub-wrapper $threads $operation $outdir $p $inpdir/$branch "${gsamples[*]}") )
      done
      
    elif [ $method == "by-branch" ]; then
      #
      # no looping necessary
      #
      outdir=$outpref.$pname/$branch
      jid+=( $(scripts-qsub-wrapper $threads $operation $outdir $p $inpdir/$branch) )

    else
      #
      # error
      #
      scripts-send2err "Error: unknown grouping method."
      exit
    fi

  done

done
       
# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



