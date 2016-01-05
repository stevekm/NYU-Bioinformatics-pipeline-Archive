#!/bin/tcsh

source ./inputs/params/params.tcsh

set aligner = bowtie2
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
set genome_index = inputs/genomes/$genome/bowtie2.index/genome
set align_params = "--local -x $genome_index"
set mapq = 30
set release = inputs/release

