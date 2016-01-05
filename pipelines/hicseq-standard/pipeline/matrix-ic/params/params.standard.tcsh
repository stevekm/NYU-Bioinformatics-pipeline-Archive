#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc               # TODO: remove all these when mirnylib is properly installed
module unload python
module load python/2.7.3

set chrom_excluded = 'chr[MY]'       # excluded chromosomes
set cutoff = 0.05

