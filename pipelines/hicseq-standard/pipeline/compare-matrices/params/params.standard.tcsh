#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'       # excluded chromosomes
set compare_params = "--min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0"           # only used if estimation was done with max-lambda=Inf

