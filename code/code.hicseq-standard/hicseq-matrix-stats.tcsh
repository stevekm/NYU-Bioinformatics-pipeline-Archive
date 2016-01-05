#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-stats.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH [OBJECTS]
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

# compute stats
set sample_paths = `echo "$branch\t$objects" | tools-key-expand | tr '\t' '/'`
Rscript ./code/hic-matrix.r stats -v -o $outdir --bin-size=$bin_size $sample_paths

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


