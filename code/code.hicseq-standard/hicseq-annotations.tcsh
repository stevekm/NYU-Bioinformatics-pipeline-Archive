#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-estimate.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH SAMPLE
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set sample = $4

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch/$sample "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# run annotations
set inpdir = $branch/$sample
set filter_branch = ../filter/results/`echo $branch | sed 's/.*\/\(filter\.\)/\1/'`
set n_reads = `cat $filter_branch/$sample/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM' | grep chr19`)
  send2err "Processing input matrix $inpdir/$mat..."
  Rscript ./code/hic-matrix.r loops -v -o $outdir -L $sample --n-reads=$n_reads $loop_params $inpdir/$mat
  send2err "Rscript ./code/hic-matrix.r loops -v -o $outdir -L $sample --n-reads=$n_reads $loop_params $inpdir/$mat"
end

# save variables
set >! $outdir/job.vars.tsv


