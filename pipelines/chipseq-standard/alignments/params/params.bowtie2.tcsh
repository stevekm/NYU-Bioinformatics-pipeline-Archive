#!/bin/tcsh

source ./inputs/params/params.tcsh

module load bowtie2
set aligner = bowtie2
set genome_index = inputs/release/../bowtie2.index/genome
set align_params = "--local -x $genome_index"
set mapq = 30
set release = inputs/release

