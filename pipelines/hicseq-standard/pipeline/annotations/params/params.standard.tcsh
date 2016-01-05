#!/bin/tcsh

source ./inputs/params/params.tcsh

set genes_bed = $genome_dir/gene.bed       # gene BED6 file for annotation of interactions
set loci_bed = ()                          # other loci of interest to be used for annotating, e.g. ChIP-seq BED files


