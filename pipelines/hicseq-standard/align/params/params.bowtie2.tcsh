#!/bin/tcsh

source ./inputs/params/params.tcsh

module load bowtie2
set aligner = bowtie2
set release = inputs/release
set genome_index = $release/../bowtie2.index/genome
set align_params = "--very-sensitive-local --local --reorder -x $genome_index"

