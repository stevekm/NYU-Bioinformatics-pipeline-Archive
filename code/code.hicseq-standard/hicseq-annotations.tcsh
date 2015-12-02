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
set sample = $4          # TODO: allow multiple samples, so that this operation can be run by-group; BUT, genomes must match!

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch/$sample "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# create table of interactions for each chromosome
set inpdir = $branch/$sample
set filter_branch = ../filter/results/`echo $branch | sed 's/.*\/\(filter\.\)/\1/'`
set n_reads = `cat $filter_branch/$sample/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv matrix.*.RData | grep chr19`)
  send2err "Processing input matrix $inpdir/$mat..."
  set outmatdir = `echo $mat | sed 's/\.tsv$//' | sed 's/\.RData$//'`
  Rscript ./code/hic-matrix.r loops -v -o $outdir/$outmatdir -L $sample --n-reads=$n_reads $loop_params $inpdir/$mat
end

# annotate
./code/hicseq-annotate-tables.tcsh $outdir/annotations "$outdir/matrix.*/loops.tsv" $genes_bed "$loci_bed"

# save variables
set >! $outdir/job.vars.tsv


