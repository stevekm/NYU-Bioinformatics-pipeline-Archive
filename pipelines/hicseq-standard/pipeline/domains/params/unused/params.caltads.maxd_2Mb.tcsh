#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = caltads
set chrom_excluded = 'chr[MY]'
set max_dist = 2000000
set win = `echo $max_dist/$bin_size | bc`
set caltads_params = "--resolution=$bin_size -w $win"

