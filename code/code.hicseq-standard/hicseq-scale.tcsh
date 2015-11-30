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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch/$sample "genome genome_dir enzyme bin_size"

# create path
scripts-create-path $outdir/

# setup
set inpdir = $branch/$sample
set n_reads = `./code/calc-matrix-reads.r $inpdir/matrix.*.tsv`
set res = `echo $bin_size/1000 | bc`
set features = $genome_dir/$enzyme.w=${res}kb.features.txt 

# run matrix scaling/normalizing operation
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM'`)
  scripts-send2err "Processing $mat [enzyme = $enzyme\; features = $features\; n_reads = $n_reads]..."
  Rscript ./code/hic-matrix.r normalize -v -o $outdir/$mat --n-reads=$n_reads --ignored-loci=$inpdir/ignored_loci.txt --features=$features $inpdir/$mat
end

# save variables
set >! $outdir/job.vars.tsv


