#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc               # this is necessary in order to take care of module conflicts in our system
module unload python
module load python/2.7.3

set chrom_excluded = 'chr[MY]'       # excluded chromosomes
set cutoff = 0.05

