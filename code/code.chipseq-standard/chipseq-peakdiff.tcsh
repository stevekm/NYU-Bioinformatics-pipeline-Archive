#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # shell settings (must be included in all scripts)

##
## USAGE: chipseq-peakdiff.tcsh OUTPUT-DIR PARAMETER-SCRIPT BRANCH OBJECTS-1 OBJECTS-2
##

# process command-line inputs
if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects1 = ($4)
set objects2 = ($5)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects1 $objects2" "genome genome_dir"

# set parameters
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- tool = $tool"
scripts-send2err "- diff params = $diff_params"
scripts-send2err "- annotation params = $annot_params"

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine input files
set xfiles = `echo $objects1 | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/matrix.tsv"}'`
set yfiles = `echo $objects2 | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/matrix.tsv"}'`

# create matrix
set mat = $outdir/matrix.tsv
scripts-multi-paste $xfiles $yfiles >! $mat
set nref = $#xfiles

# run easydiff
scripts-easydiff.r -v -o $outdir --nref=$nref $diff_params $mat

# merge overlapping diff-peaks
cat $outdir/diff.gain | sed 's/:/\t/' | sed 's/-/\t/' | scripts-sortbed | gtools-regions link --label-func max | cut -f-4 >! $outdir/diff.gain.bed
cat $outdir/diff.loss | sed 's/:/\t/' | sed 's/-/\t/' | scripts-sortbed | gtools-regions link --label-func min | cut -f-4 >! $outdir/diff.loss.bed

# annotate diff-peaks
( \
 echo "LOG-FOLD-CHANGE\tDIFF-PEAK-LOCUS\tREGION\tGENE-SYMBOL\tDISTANCE"; \
 cat $outdir/diff.gain.bed $outdir/diff.loss.bed | gtools-overlaps $annot_params | cut -f1-3,4,7 | sed 's/:[^:]\+|/\t/' | sort -u \
) >! $outdir/peakdiff.annotated.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



