#!/bin/tcsh
source ./code/code.main/custom-tcshrc   # shell settings

##
## USAGE: chipseq-heatmaps.tcsh OUTPUT-DIR PARAMETER-SCRIPT MATRICES-BRANCH
##

# process command-line inputs
if ($#argv != 3) then
  grep '^##' $0
  exit
endif

# inputs
set outdir = $1
set params = $2
set branch = $3

# create output path
scripts-create-path $outdir/

# parameters
scripts-send2err "Processing parameter file $params..."
source $params
set | scripts-send2err 

# process sample sheet
set sheet = inputs/sample-sheet.tsv
set groups_noinput = `cat $sheet | grep -v '^#' | cut -f2 | tr '|' '\n' | sort | uniq | grep -vwi input`
if ($#groups_noinput < 1) then
  set awk_filter = '$0'                        # include inputs only if no more than one group (besides inputs) is available
else
  set awk_filter = 'tolower($2)!~/input/'      # remove all inputs
  scripts-send2err "Warning: Excluding all Inputs from heatmaps."
endif

# create dataset file [TODO: this will fail if you have multiple group names per sample!]
cat $sheet | grep -v '^#' | awk $awk_filter | cut -f1,2 | awk -v pref=$branch '{print $2":"$1"\t"pref"/"$1"/matrix.tsv"}' >! $outdir/dataset.tsv

# create heatmaps
scripts-heatclustering.r -v -o $outdir $heatmap_params $outdir/dataset.tsv 

# cleanup
rm -f $outdir/*.RData


