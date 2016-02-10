#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'       # excluded chromosomes

set max_dist = `echo 10000000/$bin_size | bc`      # number of bins (max distance = 10Mb)

set compare_params = "--max-dist=$max_dist --n-dist=1 --min-lambda=0.0 --max-lambda=1.0 --n-lambda=6 --gamma=0"           # only used if estimation was done with max-lambda=Inf

