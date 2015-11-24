#!/bin/tcsh

# load basic tools
module unload samtools
module unload java
module load samtools/1.2.1
module load bedtools/2.22.0
module load java/1.7
module load picard-tools

set sheet = inputs/sample-sheet.tsv

