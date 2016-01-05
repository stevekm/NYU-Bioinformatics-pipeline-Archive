#!/bin/tcsh

source ./inputs/params/params.tcsh

set peakdist = 500      # max distances between peaks that will be merged
set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"

