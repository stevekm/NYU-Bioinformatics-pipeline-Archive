#!/bin/bash

##
## USAGE: scripts-master-loop.sh THREADS SCRIPT-WRAPPER OUTDIR-PREFIX PARAM-SCRIPTS INPUT-DIR
##

if (( $# != 5 )); then
  grep '^##' $0
  exit
fi

threads=$1
wrapper=$2
outpref=$3
params=($4)
inpdir=$5
groups=$6

# start
jid=()
# loop over all parameter files
for p in ${params[*]}; do
  pname=$(echo $p | sed 's/.*\///' | sed 's/^params\.//' | cut -d'.' -f1)

  # loop over all branches of input tree
  for leaf in $(cd $inpdir; find . -name 'job.sh' | sed 's/\/job.sh$//' | sed 's/^.\///'); do
    outdir=$outpref.$pname/$leaf
    jid+=( $(scripts-qsub-wrapper $threads $wrapper $outdir $p $inpdir/$leaf) )
  done

done
       
# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



