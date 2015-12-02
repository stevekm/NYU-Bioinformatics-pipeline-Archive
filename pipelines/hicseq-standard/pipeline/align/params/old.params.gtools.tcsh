#!/bin/tcsh

source ./inputs/params/params.tcsh

module load bowtie2
set aligner = gtools
set genome = `cat $sheet | awk -v s=$sample '$1==s' | cut -f5`
set genome_index = inputs/genomes/$genome/bowtie2.index/genome
set align_params = "--min-len 25 --len-diff 5"

