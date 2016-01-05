#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-estimate.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)          # TODO: allow multiple objects (i.e. samples), so that this operation can be run by-group; BUT, genomes must match!

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# annotate
set inpdir = $branch/$object
./code/hicseq-annotate-tables.tcsh $outdir "$inpdir/matrix.*/loops.tsv" $genes_bed "$loci_bed"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


