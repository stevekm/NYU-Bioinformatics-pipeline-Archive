#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-diff-domains.tcsh OUTDIR PARAMS BRANCH OBJECT1 OBJECT2
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
    Rscript ./code/hic-matrix.r domain-diff -v -o $outdir/$chr --row-labels $diff_domains_params $branch/$object1/$f $branch/$object2/$f             # TODO: use scripts-qsub-run to assign memory
  endif
end

# Collect diff-domains from all chromosomes
set K = `ls -1 $outdir/*/table.k=*.tsv | sed 's/.*\/table.k=//' | cut -d'.' -f1 | sort -u`
foreach k ($K)
  set t = table.k=$k.tsv
  cat $outdir/*/$t | head -1 >! $outdir/$t
  foreach tt ($outdir/*/$t)
    cat $tt | sed '1d' >> $outdir/$t
  end
  cat $outdir/$t | sed '1d' | cut -f-2 | awk '$2==1' | cut -f1 | tr ':-' '\t' | gtools-regions shiftp -5p -1 -3p 0 >! $outdir/boundary_gain.k=$k.bed
  cat $outdir/$t | sed '1d' | cut -f-2 | awk '$2==-1' | cut -f1 | tr ':-' '\t' | gtools-regions shiftp -5p -1 -3p 0 >! $outdir/boundary_loss.k=$k.bed  
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




