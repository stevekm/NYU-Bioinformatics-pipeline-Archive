#!/bin/tcsh

##
## USAGE: chipseq-pca.tcsh OUTPUT-DIR PARAMETER-SCRIPT MATRICES-BRANCH
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 3) then
  grep '^##' $0
  exit
endif

# inputs
set out = $1
set params = $2
set branch = $3

# create output path
scripts-create-path $out/

# create record of input branch(es) [TODO: need to figure out a way to do a relative symlink]
echo $branch > $out/obj.branch

# parameters
scripts-send2err "Processing parameter file $params..."
source $params

# process sample sheet
set sheet = inputs/sample-sheet.tsv
set groups_noinput = `cat $sheet | grep -v '^#' | cut -f2 | tr '|' '\n' | sort | uniq | grep -vwi input`
if ($#groups_noinput < 1) then
  set awk_filter = '$0'                        # include inputs only if no more than one group (besides inputs) is available
else
  set awk_filter = 'tolower($2)!~/input/'      # remove all inputs from PCA
  scripts-send2err "Warning: Excluding all Inputs from PCA."
endif

# determine labels and matrices [# TODO: the code below will not work for multi-group assignments]
set labels = `cat $sheet | grep -v '^#' | awk $awk_filter | cut -f1,2 | cut -d'|' -f1 | awk '{print $2":"$1}'`
set matrices = `cat $sheet | grep -v '^#' | awk $awk_filter | cut -f1 | awk -v pref=$branch '{print pref"/"$1"/matrix.tsv"}'`

# run PCA
echo $labels | tr ' ' '\n' >! $out/labels.tsv
scripts-multi-paste $matrices | grep -v '^chr[MYX]' >! $out/matrix.tsv     # TODO: exclude chrM/X/Y?
scripts-perform-pca.r -v -o $out -L $out/labels.tsv --show-text --use-short-names $out/matrix.tsv

# cleanup
rm -f $out/matrix.tsv




