#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hic-filter.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-BRANCH SAMPLE
##

# process command-line inputs
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

# indentify genome directory
set genome = `cat $sheet | awk -v s=$sample '$1==s' | cut -f5`
set enzyme = `cat $sheet | awk -v s=$sample '$1==s' | cut -f6`
set genome_dir = inputs/genomes/$genome

# create path
scripts-create-path $outdir/

# filter
samtools view $branch/$sample/alignments.bam | gtools_hic filter -v -E $genome_dir/$enzyme.fragments.bed $filter_params | gzip >! $outdir/filtered.reg.gz

# save variables
set >! $outdir/job.vars.tsv





