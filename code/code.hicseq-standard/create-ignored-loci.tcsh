#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: create-ignored-loci.tcsh GENOME-DIR BIN-SIZE
##

if ($#argv != 2) then
  grep '^##' $0
  exit
endif

set genome_dir = $1
set bin_size = $2

cat $genome_dir/genome.bed | gtools-regions win -s $bin_size -d $bin_size | gtools-overlaps subset -i $genome_dir/centrotelo.bed | gtools-regions shiftp -5p 1 -3p 0 | sed 's/\t/:/' | sed 's/\t/-/' | cut -f1


