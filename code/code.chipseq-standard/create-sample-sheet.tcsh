#!/bin/tcsh

##
## USAGE: create-sample-sheet
##
## FUNCTION: create sample sheet automatically from input files in fastq-or-alignments directory
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# check for duplicate sample names
set dups = `(cd fastq-or-alignments; ls -1 *.fastq *.fastq.gz *.bam) | sed 's/\.gz$//' | sed 's/\.fastq$//' | sed 's/\.bam$//' | sort | uniq -d | wc -l` 
if ($dups > 0) then
  scripts-send2err "Error: duplicate sample names, check your fastq and bam files!"
  exit
endif

set sheet = sample-sheet.tsv
echo "#SAMPLE-NAME\tGROUPS\tCONTROL-SAMPLE-NAME\tFILE-NAME" >! $sheet
set files = `cd fastq-or-alignments; ls -1 *.fastq *.fastq.gz *.bam`
foreach f ($files)
  set sample = `echo $f | sed 's/\.gz$//' | sed 's/\.fastq$//' | sed 's/\.bam$//'`
  set group = `echo $sample | cut -d'-' -f-3`
  ( \
    echo -n $sample'\t' ; \
    echo -n $group'\t' ; \
    echo -n 'n/a\t' ; \
    echo $f \
  ) >> $sheet
end

cat $sheet

echo "tissue/cell-line: "
grep -v '^#' $sheet | cut -f1 | cut -d'-' -f1 | sort | uniq -c
echo "treatment: "
grep -v '^#' $sheet | cut -f1 | cut -d'-' -f2 | sort | uniq -c
echo "ChIP: "
grep -v '^#' $sheet | cut -f1 | cut -d'-' -f3 | sort | uniq -c

