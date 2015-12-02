#!/bin/tcsh

source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-tracks.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH SAMPLES
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set samples = ($4) 

# sample paths
set sample_paths = `echo "$branch\t$samples" | key_expand | tr '\t' '/'`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars "$sample_paths" "genome genome_dir enzyme"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# Plot barplots
set sample_reads = `echo $sample_paths | tr ' ' '\n' | sed 's/$/\/filtered.reg.gz/'`
./code/hicseq-tracks-washu.tcsh $outdir "$sample_reads" $genome_dir/genome.bed $bin_size

# save variables
set >! $outdir/job.vars.tsv

