#!/bin/tcsh

source ./inputs/params/params.tcsh

# Load required modules
module unload gcc
module load matlab/R2013a

set tool = di
set domaincallpath = "/ifs/home/cl3011/ROTATION_3/Resources/Software/domaincall_software"
set chrom_excluded = 'chr[MY]'
set window_size = 2000000 #2MB upstream and 2MB downstream 
set min = 2
set prob = 0.99
set faipath = $genome_dir/bowtie2.index/genome.fa.fai



