#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'       # excluded chromosomes

set loop_params = "--bin-size=$bin_size --lambda-id=6 --rpk2b-cutoff=1.0 --loop-cutoff=4.0 --min-distance=40000"        # parameters for identifying significant interactions


