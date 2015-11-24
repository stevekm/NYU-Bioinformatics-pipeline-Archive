#!/bin/tcsh

##
## USAGE: create-sample-sheet.tcsh GENOME={hg19,mm10}
##
## FUNCTION: create sample sheet automatically from input files in fastq directory
##

# shell settings
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 1) then
  grep '^##' $0
  exit
endif

# inputs/outputs
set genome = $1
set inpdir = fastq
set sheet = sample-sheet.tsv

# check for duplicate sample names
set dups = `(cd $inpdir; ls -1 *.fastq *.fastq.gz) | sed 's/\.gz$//' | sed 's/\.fastq$//' | sort | uniq -d | wc -l` 
if ($dups > 0) then
  scripts-send2err "Error: duplicate sample names, check your fastq and bam files!"
  exit
endif

set files = `scripts-create-temp`
set enzymes = `scripts-create-temp`
( cd $inpdir ; ls -1 *.fastq.gz ) | cols 0 0 | sed 's/_R[12].fastq.gz /\t/' | sort | mergeuniq -merge | sed 's/ *$//' | tr ' ' '\t' >! $files
cut -f1 $files | sed 's/.*HindIII.*/HindIII/' | sed 's/.*NcoI.*/NcoI/' >! $enzymes
echo "#SAMPLE-NAME\tGROUP-NAMES(S)\tFASTQ-R1-FILES\tFASTQ-R2-FILES\tGENOME\tENZYME" >! $sheet
cut -f1 $files | cut -d'-' -f-2 | paste - $files | cols -t 1 0 2 3 | sed 's/$/\t'$genome'/' | paste - $enzymes >> $sheet
rm -f $files $enzymes

echo "Your sample sheet has been created! Here is how it looks:"
echo
cat $sheet
echo

echo "Diagnostics: "
echo "Field #1: "
grep -v '^#' $sheet | cut -f1 | cut -d'-' -f1 | sort | uniq -c
echo "Field #2: "
grep -v '^#' $sheet | cut -f1 | cut -d'-' -f2 | sort | uniq -c
echo "Field #3: "
grep -v '^#' $sheet | cut -f1 | cut -d'-' -f3 | sort | uniq -c

