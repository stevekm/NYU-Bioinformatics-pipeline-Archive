#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-boundary-scores.tcsh OUTDIR PARAMS MATRIX-BRANCH OBJECTS
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
  scripts-send2err "Error: hicseq-boundary-scores.tcsh requires exactly one input object."
  exit 1
endif

set object1 = $objects[1]
foreach f (`cd $branch/$object1; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`)
  set chr = `echo $f | cut -d'.' -f2`
  scripts-send2err "Processing matrix $f..."
  Rscript ./code/hic-matrix.r domains -v -o $outdir/$chr --row-labels $boundary_scores_params $branch/$object1/$f
end

# Collect boundary scores from all chromosomes
set pref = all_scores
set K = `ls -1 $outdir/*/$pref.k=*.tsv | sed "s/.*\/$pref.k=//" | cut -d'.' -f1 | sort -u`
foreach k ($K)
  set t = $pref.k=$k.tsv
  cat $outdir/*/$t | head -1 >! $outdir/$t
  foreach tt ($outdir/*/$t)
    cat $tt | sed '1d' >> $outdir/$t
  end
end

# remove unused files
rm -rf $outdir/*/

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




