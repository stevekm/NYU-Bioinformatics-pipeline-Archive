#!/bin/tcsh

source ./inputs/params/params.tcsh

module load bowtie2
set aligner = bowtie2
set genome = `cat $sheet | awk -v s=$sample '$1==s' | cut -f5`
set genome_index = inputs/genomes/$genome/bowtie2.index/genome
set align_params = "--very-sensitive-local --local --reorder"

