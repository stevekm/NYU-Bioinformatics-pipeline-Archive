#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains-hicmatrix.tcsh OUTDIR PARAMS BRANCH OBJECTS
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
  scripts-send2err "Error: hicseq-domains-hicmatrix.tcsh requires exactly one input object."
  exit 1
endif

set object1 = $objects[1]
foreach f (`cd $branch/$object1; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`)
  set chr = `echo $f | cut -d'.' -f2`
  scripts-send2err "Processing matrix $f..."
  Rscript ./code/hic-matrix.r domains -v -o $outdir/$chr --row-labels $hicmatrix_params $branch/$object1/$f
end

# Collect domains from all chromosomes
set pref = boundaries
set K = `ls -1 $outdir/*/$pref.k=*.tsv | sed "s/.*\/$pref.k=//" | cut -d'.' -f1 | sort -u`
foreach k ($K)
  set t = $pref.k=$k.tsv
  cat $outdir/*/$pref.k=$k.tsv | sed 's/:/\t/' | sed 's/-/\t/' | gtools-regions shift -start -1 -stop 0 | gtools-regions inv -g $genome_dir/genome.bed | cut -f-3 >! $outdir/domains.k=$k.bed
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




