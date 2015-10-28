#!/bin/tcsh

##
## USAGE: hicseq-align.tcsh OUTPUT-DIR PARAM-SCRIPT FASTQ1 FASTQ2
##

# shell settings
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set fastq1 = ($3)
set fastq2 = ($4)

# set parameters
source $params
if (! $?NSLOTS) then
  set threads = 8
else
  set threads = $NSLOTS
endif

# create path
scripts-create-path $out

# align
if ($aligner == 'gtools') then                    ## Aligner = gtools

  if ( ($#fastq1 > 1) || ($#fastq2 > 1)) then
    script-send2err "Error: gtools_hic align does not allow multiple read pair files."
    exit
  endif
  gtools_hic align -v --work-dir $out/tmp -p $threads $align_params $fastq1 $fastq2 | samtools view -T $release/../bowtie2.index/genome.fa -b1 - >! $out/alignments.bam
  rm -rf $out/tmp


else if ($aligner == 'bowtie2') then              ## Aligner = bowtie2

  bowtie2 -p $threads $align_params --sam-nohead -U `echo $fastq1 | tr ' ' ','` -S $out/alignments.R1.sam
  bowtie2 -p $threads $align_params --sam-nohead -U `echo $fastq2 | tr ' ' ','` -S $out/alignments.R2.sam
  paste -d'\n' $out/alignments.R1.sam $out/alignments.R2.sam | samtools view -T $release/../bowtie2.index/genome.fa -b1 - >! $out/alignments.bam
  rm -f $out/alignments.R1.sam $out/alignments.R2.sam


else
  script-send2err "Error: unknown aligner $aligner."
  exit
endif





