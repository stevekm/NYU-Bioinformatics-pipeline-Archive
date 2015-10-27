#!/bin/tcsh

##
## USAGE: create-sample-sheet
##
## FUNCTION: create sample sheet automatically from input files in fastq-or-alignments directory
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# set path
set path = (./code/code $path)

# inputs/outputs
set inpdir = fastq
set sheet = sample-sheet.tsv

# check for duplicate sample names
set dups = `(cd $inpdir; ls -1 *.fastq *.fastq.gz *.bam) | sed 's/\.gz$//' | sed 's/\.fastq$//' | sed 's/\.bam$//' | sort | uniq -d | wc -l` 
if ($dups > 0) then
  scripts-send2err "Error: duplicate sample names, check your fastq and bam files!"
  exit
endif

echo "#SAMPLE-NAME\tFASTQ-R1-FILES\tFASTQ-R2-FILES\tGROUPS" >! $sheet
( cd $inpdir ; ls -1 *.fastq.gz ) | cols 0 0 | sed 's/_R[12].fastq.gz / /' | sort | mergeuniq -merge | sed 's/ *$//' | tr ' ' '\t' >> $sheet

cat $sheet

echo "Field #1: "
grep -v '^#' $sheet | cut -f2 | cut -d'-' -f1 | sort | uniq -c
echo "Field #2: "
grep -v '^#' $sheet | cut -f2 | cut -d'-' -f2 | sort | uniq -c
echo "Field #3: "
grep -v '^#' $sheet | cut -f2 | cut -d'-' -f3 | sort | uniq -c

