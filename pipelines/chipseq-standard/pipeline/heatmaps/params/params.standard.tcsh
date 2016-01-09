#!/bin/tcsh

source ./inputs/params/params.tcsh

set group = 'group'             # sample sheet column name used to assign group colors
set heatmap_params = "--use-short-names --normalize=none --log2=true --n-best=4000 --diff=0.75 --nclust=10"


