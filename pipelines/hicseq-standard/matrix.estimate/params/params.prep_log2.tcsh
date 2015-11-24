#!/bin/tcsh

source ./inputs/params/params.tcsh

set prep = log2
set max_lambda = 1.0
set n_lambda = 6
set log2_lambda = #'--log2-lambda'
if ($bin_size >= 1000000) then
  set algorithm = fused2D_flsa
  set zone_size = 0
else
  set algorithm = fused2Dzone_flsa
  set zone_size = `echo 4000000/$bin_size | bc`        # number of bins
endif
set skip_dist = 0
set split_size = 100
set gamma = (0.00)

set hic_params = "--flsa-verbose --row-labels --skip-distance=$skip_dist --preprocess=$prep --min-lambda=0 --max-lambda=$max_lambda --n-lambda=$n_lambda $log2_lambda --algorithm=$algorithm --zone-size=$zone_size --split-check-size=$split_size"



