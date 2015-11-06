#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # custom shell environment

##
## USAGE: validate-sample-sheet SAMPLE-SHEET-TSV
##

if ($#argv != 1) then
  grep '^##' $0
  exit 1
endif

set sheet = $1

# check for SPACES
if (`cat $sheet | grep -v '^#' | grep ' ' | wc -l` > 0) then
  scripts-send2err 'Error: spaces are not allowed in sample sheet.'
  exit 1
endif

# duplicate sample names (column #1)
if (`cat $sheet | grep -v '^#' | cut -f1 | sort | uniq -d | wc -l` > 0) then
  scripts-send2err 'Error: duplicate sample names are not allowed in column 1 of the sample sheet.'
  exit 1
endif

# control sample names (column #3) not found in sample names (column #1)
set t = `scripts-create-temp`
cat $sheet | grep -v '^#' | cut -f3 | grep -v '^n/a$' | sort -u >! $t
if (`cat $sheet | grep -v '^#' | cut -f1 | sort | join -v1 $t - | wc -l` > 0) then
  scripts-send2err 'Error: control sample names in column 3 not found in sample names in column 1.'
  cat $sheet | grep -v '^#' | cut -f1 | sort | join -v1 $t - | scripts-send2err
  rm -f $t
  exit 1
endif
rm -f $t

# file names should exist in fastq-or-alignments directory
set files = `cat $sheet | grep -v '^#' | cut -f4 | tr ',' '\n' | sort -u`
foreach f ($files)
  if (! -e inputs/fastq-or-alignments/$f) then
    scripts-send2err "Error: file \'$f\' not found in inputs/fastq-or-alignments directory."
    exit 1
  endif
end



