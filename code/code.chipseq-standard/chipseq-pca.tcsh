#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-pca.tcsh OUTPUT-DIR PARAM-SCRIPT MATRICES-BRANCH [OBJECTS]
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

# determine labels and matrices
set labels = `cat $sheet | sed '1d' | cut -f1,$group_column | awk $awk_filter | awk '{print $2":"$1}'`
set matrices = `cat $sheet | sed '1d' | cut -f1,$group_column | awk $awk_filter | cut -f1 | awk -v pref=$branch '{print pref"/"$1"/matrix.tsv"}'`

# run PCA
echo $labels | tr ' ' '\n' >! $outdir/labels.tsv
scripts-multi-paste $matrices | grep -vwE "$chrom_excluded" >! $outdir/matrix.tsv
scripts-perform-pca.r -v -o $outdir -L $outdir/labels.tsv --show-text --use-short-names $outdir/matrix.tsv

# cleanup
rm -f $outdir/matrix.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




