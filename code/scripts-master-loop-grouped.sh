#!/bin/bash

##
## USAGE: scripts-master-loop-grouped.sh THREADS SCRIPT-WRAPPER OUTDIR-PREFIX PARAM-SCRIPTS INPUT-DIR GROUPS
##

if (( $# != 6 )); then
  grep '^##' $0
  exit
fi

threads=$1
wrapper=$2
outpref=$3
params=($4)
inpdir=$5
groups=($6)

# set path
PATH=./code/code:$PATH

# start
sheet=inputs/sample-sheet.tsv
jid=()
# loop over all parameter files
for p in ${params[*]}; do
  pname=$(echo $p | sed 's/.*\///' | sed 's/^params\.//' | cut -d'.' -f1)

  # loop over all branches of input tree (up to one level before the leaves)
  for leaf2 in $(cd $inpdir; find . -name 'job.sh' | sed 's/\/job.sh$//' | sed 's/^.\///' | sed 's/\/[^/]\+$//' | sort -u); do
    # loop over all sample groups
    for grp in ${groups[*]}; do
      outdir=$outpref.$pname/$leaf2/$grp
      leaves=( $(cat $sheet | awk -v grp=$grp '$4==grp' | cut -f1 | awk -v d=$inpdir/$leaf2 '{print d"/"$0}') )      # TODO: need to modify this for samples with multiple group assignments
      jid+=( $(scripts-qsub-wrapper $threads ./code/hicseq-matrix.tcsh $outdir $params "${leaves[*]}") )
    done
  done

done
       
# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



