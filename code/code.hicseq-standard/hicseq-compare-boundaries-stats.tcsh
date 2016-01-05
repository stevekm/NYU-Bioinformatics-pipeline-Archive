#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compare-boundaries-stats.tcsh OUTDIR PARAM-SCRIPT DOMAIN-BRANCH [OBJECTS]
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

# Check for number of objects
if ($#objects < 2) then
	scripts-send2err "Error: more than one input objects are required."
	exit 1
endif 

# Collect results and create plots
set files = `echo "$branch\t$objects" | tools-key-expand | tr '\t' '/' | sed 's/$/\/comparison.tsv/'`
head -1 $files[1] >! $outdir/comparisons.tsv
foreach x ($files)
  cat $x | scripts-skipn 1 >> $outdir/comparisons.tsv
end
Rscript ./code/correlogram-boundaries.r $outdir/comparisons.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




