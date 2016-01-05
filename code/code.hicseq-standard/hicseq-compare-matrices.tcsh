#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compare-matrices.tcsh OUTDIR PARAMS BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = $4
set object2 = $5

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$object1 $object2" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

foreach f (`cd $branch/$object1; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`)
  set chr = `echo $f | cut -d'.' -f2`
  scripts-send2err "Processing matrix $f..."
  if ((-e $branch/$object1/$f) && (-e $branch/$object2/$f)) then
    Rscript ./code/hic-matrix.r compare -v -o $outdir/$chr $compare_params $branch/$object1/$f $branch/$object2/$f 
  endif
end

# Collect all the pearson and spearman coefficients along with sample and lambda info
set comp = `basename $outdir`
set methods = "pearson spearman"
foreach method ($methods)
	cat $outdir/*.cor.$method.tsv | grep -v lambda | sed "s/^/$object1	$object2	${comp}	${method}	/" >! $outdir/cor.$method.tsv
end 

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




