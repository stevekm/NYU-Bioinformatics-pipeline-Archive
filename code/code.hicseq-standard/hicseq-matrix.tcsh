#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix.tcsh OUTPUT-DIR PARAM-SCRIPT FILTER-BRANCH SAMPLES
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set samples = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch/$samples[1] "genome genome_dir enzyme"                      # TODO: check consistency of these variables across samples

# run parameter script
scripts-send2err "Setting parameters..."
source $params

# create path
scripts-create-path $outdir/

# setup
cat $genome_dir/genome.bed | genomic_regions win -s $bin_size -d $bin_size | genomic_overlaps subset -i $genome_dir/centrotelo.bed | genomic_regions shiftp -5p 1 -3p 0 | sed 's/\t/:/' | sed 's/\t/-/' >! $outdir/ignored_loci.txt

# generate matrix
set filtered_reads = `echo $samples | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/filtered.reg.gz"}'`
scripts-smartcat $filtered_reads | gtools_hic matrix -v -p $outdir/matrix. $matrix_params

# save variables
set >! $outdir/job.vars.tsv

