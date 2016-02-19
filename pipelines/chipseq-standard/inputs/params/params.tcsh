#!/bin/tcsh

# load all tools
module unload gcc
module unload samtools
module unload java
module unload r

module load samtools/1.2.1
module load bedtools/2.22.0
module load java/1.7
module load picard-tools
module load bowtie2
#module load macs/2.0.10.20131216
module load r/3.2.0

# set sample sheet
set sheet = inputs/sample-sheet.tsv

