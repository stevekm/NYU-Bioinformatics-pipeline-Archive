#!/bin/tcsh

source ./inputs/params/params.tcsh

# TODO: this should be run for all lambdas, need to modify hic-matrix.r loops operation
set loop_params = "--bin-size=$bin_size --lambda-id=6 --rpk2b-cutoff=1.0 --loop-cutoff=4.0 --min-distance=40000"        # parameters for identifying significant interactions
set genes_bed = $genome_dir/gene.bed       # gene BED6 file for annotation of interactions
set loci_bed = ()                          # other loci of interest to be used for annotating, e.g. ChIP-seq BED files


