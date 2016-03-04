#!/bin/tcsh

# load basic tools
module unload samtools
module unload java
module unload gcc
module unload python
module load samtools/1.2.1
module load bedtools/2.22.0
module load java/1.7
module load picard-tools
module load r/3.0.2

# load tools required for each step of the pipeline (this can be overriden in local param scripts)
module load bowtie2/2.2.6
module load armatus/2014-05-19
module load caltads/0.1.0
module load ghmm/0.9

# sample sheet file
set sheet = inputs/sample-sheet.tsv

