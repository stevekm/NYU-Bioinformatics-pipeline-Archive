#!/bin/tcsh

source ./inputs/params/params.tcsh

set group = 'group'                            # sample sheet column name used to assign group colors
set include_input = 'false'                    # include ChIP inputs
set chrom_excluded = 'chr[MXY]'                # excluded chromosomes
set heatmap_params = "--use-short-names --normalize=mean --log2=true --n-best=4000 --diff=0.75 --nclust=10"

