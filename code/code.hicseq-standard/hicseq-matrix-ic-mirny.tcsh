#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-ic-mirny.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S) MATRIX
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)          # TODO: allow multiple samples, so that this operation can be run by-group
set mat = $5

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# convert to IC input format
scripts-send2err "Running iterative correction..."
set inpdir = $branch/$object
set x = $outdir/$mat.inp
set y = $outdir/$mat.out
cat $inpdir/$mat | scripts-skipn 1 | cut -f2- >! $x
  
# run IC
python ./code/hicseq-matrix-ic-mirny.py $x $cutoff $y

# convert IC output to standard output matrix
scripts-send2err "Converting to standard matrix format..."
cat $inpdir/$mat | head -1 >! $outdir/$mat
cat $inpdir/$mat | scripts-skipn 1 | cut -f1 | paste - $y | tr -s ' ' '\t' >> $outdir/$mat

# cleanup
rm -f $x $y

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


