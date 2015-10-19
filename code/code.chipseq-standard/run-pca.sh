#!/bin/bash

##
## USAGE: run.pca
##

if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# set path
PATH=./code/code:$PATH

# process sample sheet
sheet=inputs/sample-sheet.tsv                                                                         # TODO: the code below will not work for multi-group assignments
groups_noinput=( $(cat $sheet | grep -v '^#' | cut -f4 | tr '|' '\n' | sort | uniq | grep -vwi input) )
if [ ${#groups_noinput[@]} -le 1 ]; then
  awk_filter='$0'                        # include inputs only if no more than one group (besides inputs) is available
else
  awk_filter='tolower($2)!~/input/'      # remove all inputs from PCA
  scripts-send2err "Warning: Excluding all Inputs from PCA."
fi

# run PCA
scripts-send2err "=== Running PCA ============="
scripts-create-path results/
threads=1
inpdir=matrices/results
D=( $(cd $inpdir; find . -name 'matrix.*.tsv' | sed 's/\/matrix.*.tsv$//' | sort -u | grep 'nbins=1/' | sed 's/^\.\///') )
jid=()
for d in ${D[*]}; do
  outdir=results/`echo $d | sed "s/matrices\./pca./"`
  scripts-send2err "- out = $outdir"
  labels=( $(cat $sheet | grep -v '^#' | awk $awk_filter | cut -f2,4 | awk '{print $2":"$1}') )
  files=( $(cat $sheet | grep -v '^#' | awk $awk_filter | cut -f2 | awk -v pref=$inpdir/$d/ '{print pref"matrix."$1".tsv"}') )
  jid+=( $(scripts-qsub-wrapper $threads ./code/matrix-pca $outdir "${labels[*]}" "${files[*]}") )
done

# done
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "$jid"
scripts-send2err "Done."





