#!/bin/tcsh

# load all tools
module unload samtools
module unload java
module load samtools/1.2.1
module load bedtools/2.22.0
module load java/1.7
module load picard-tools
module load bowtie2
#module load macs/2.0.10.20131216

# set sample sheet
set sheet = inputs/sample-sheet.tsv

