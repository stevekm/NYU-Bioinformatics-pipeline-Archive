#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hic-filter.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# run parameter script
source $params

# indentify genome directory
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
set enzyme = `./code/read-sample-sheet.tcsh $sheet $object enzyme`
set genome_dir = inputs/genomes/$genome

# create path
scripts-create-path $outdir/

# filter
samtools view $branch/$object/alignments.bam | gtools-hic filter -v -E $genome_dir/$enzyme.fragments.bed --stats $outdir/stats.tsv $filter_params | gzip >! $outdir/filtered.reg.gz

# save variables
set >! $outdir/job.vars.tsv





