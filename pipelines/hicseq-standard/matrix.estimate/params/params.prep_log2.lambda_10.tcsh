#!/bin/tcsh

source ./inputs/params/params.tcsh

set res = 1000
set prep = log2
set max_lambda = 10
set n_lambda = 5
set log2_lambda = '--log2-lambda'
set algorithm = fused2D_flsa
set zone_size = `echo 2000/$res | bc`        # number of bins
set skip_dist = 0
set split_size = 100
set gamma = (0.00)

set hic_params = "--flsa-verbose --row-labels --skip-distance=$skip_dist --preprocess=$prep --min-lambda=0 --max-lambda=$max_lambda --n-lambda=$n_lambda $log2_lambda --algorithm=$algorithm --zone-size=$zone_size --split-check-size=$split_size"



