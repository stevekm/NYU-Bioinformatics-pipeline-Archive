#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-annotations-stats.tcsh OUTPUT-DIR PARAM-SCRIPT ANNOTATIONS-BRANCH OBJECTS
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

if ($#objects != 1) then
  scripts-send2err "Error: hicseq-annotations-stats.tcsh requires exactly one input object."
  exit 1
endif

# compute stats
Rscript ./code/hicseq-annotations-enrichments.r $outdir $branch/$objects[1]/table.annotated.tsv $nbest

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


