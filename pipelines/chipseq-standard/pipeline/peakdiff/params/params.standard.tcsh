#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = "easydiff"
set diff_params = "--normalize=none --scale=log2 --method=mean --fc-cutoff=2.0 --val-cutoff=0"
set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"


