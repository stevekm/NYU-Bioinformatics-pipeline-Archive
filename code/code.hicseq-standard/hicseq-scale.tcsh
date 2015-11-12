#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-scale.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH SAMPLE
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set sample = $4

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# setup
set inpdir = $branch/$sample
set n_reads = `cat $inpdir/matrix.*.tsv | cut -f2- | grep -v '^chr' | tr '\t' '\n' | vectors sum -n 0 | matrix csum -l`                 # TODO: this is not correct because the off-diagonal elements should not be counted twice!
set bin_size = `cat $inpdir/ignored_loci.txt | head -1 | sed 's/:/ + /' | sed 's/-/ /' | sed 's/^/_\t/' | genomic_regions n | cut -f2`
set res = `echo $bin_size/1000 | bc`
set enzyme = HindIII
if (`echo $sample | grep -c 'NcoI'` == 1) then
  set enzyme = NcoI
endif
set features = inputs/release/../DNA/features/$enzyme.w=${res}kb.features.txt 

# run matrix scaling/normalizing operation
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM'`)
  scripts-send2err "Processing $mat [enzyme = $enzyme\; features = $features\; n_reads = $n_reads]..."
  ./code/hic_matrix.r normalize -v -o $outdir/$mat --n-reads=$n_reads --ignored-loci=$inpdir/ignored_loci.txt --features=$features $inpdir/$mat
end


