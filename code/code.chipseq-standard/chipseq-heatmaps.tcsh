#!/bin/tcsh
source ./code/code.main/custom-tcshrc   # shell settings

##
## USAGE: chipseq-heatmaps.tcsh OUTPUT-DIR PARAMETER-SCRIPT MATRICES-BRANCH [OBJECTS]
##

# process command-line inputs
if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

# inputs
set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if samples is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# process sample sheet
set group_column = `cat $sheet | head -1 | tr '\t' '\n' | grep -n '^group$' | cut -d':' -f1`
set groups_noinput = `cat $sheet | sed '1d' | cut -f$group_column | sort -u | grep -vwi input`
if ($#groups_noinput == 0) then
  set awk_filter = '$0'                        # include inputs only if no more than one group (besides inputs) is available
else
  set awk_filter = 'tolower($2)!~/input/'      # remove all inputs from PCA
  scripts-send2err "Warning: Excluding all Inputs from PCA."
endif

# create dataset file
cat $sheet | sed '1d' | cut -f1,$group_column | awk $awk_filter | awk -v pref=$branch '{print $2":"$1"\t"pref"/"$1"/matrix.tsv"}' >! $outdir/dataset.tsv

# create heatmaps
scripts-heatclustering.r -v -o $outdir $heatmap_params $outdir/dataset.tsv 

# cleanup
rm -f $outdir/*.RData

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


