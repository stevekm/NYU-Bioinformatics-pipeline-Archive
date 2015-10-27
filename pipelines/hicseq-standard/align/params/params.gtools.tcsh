#!/bin/tcsh

source ./inputs/params/params.tcsh

module load bowtie2
set aligner = gtools
set release = inputs/release
set genome_index = $release/../bowtie2.index/genome
set align_params = "--min-len 25 --len-diff 5 --bowtie-index $genome_index"

