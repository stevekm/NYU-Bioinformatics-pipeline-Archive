#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compare-matrices-stats.tcsh OUTDIR PARAM-SCRIPT COMPARE-MATRICES-BRANCH [OBJECTS]
##

if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

scripts-send2err "Generating correlograms..."
set methods = "pearson spearman"
foreach method ($methods)
  scripts-create-path $outdir/$method/ 
  set files = `echo "$branch\t$objects" | tools-key-expand | tr '\t' '/' | sed 's/$/\/'"cor.$method.tsv/"`
  cat $files >! $outdir/$method/cor.$method.tsv
  Rscript ./code/correlogram-matrices.r $outdir/$method/cor.$method.tsv
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




