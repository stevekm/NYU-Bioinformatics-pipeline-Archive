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
if (`cat $sheet | sed '1d' | grep ' ' | wc -l` > 0) then
  scripts-send2err 'Error: spaces are not allowed in sample sheet.'
  exit 1
endif

# duplicate sample names (column #1)
if (`cat $sheet | sed '1d' | cut -f1 | sort | uniq -d | wc -l` > 0) then
  scripts-send2err 'Error: duplicate sample names are not allowed in column 1 of the sample sheet.'
  exit 1
endif

# control sample names (column #2) not found in sample names (column #1)
set control_col = `cat $sheet | head -1 | tr '\t' '\n' | grep -n '^control$' | cut -d':' -f1`
set t = `scripts-create-temp`
cat $sheet | sed '1d' | cut -f$control_col | grep -v '^NA$' | sort -u >! $t
if (`cat $sheet | sed '1d' | cut -f1 | sort | join -v1 $t - | sort -u | wc -l` > 0) then
  scripts-send2err 'Error: control sample names in column $control_col not found in sample names in column 1.'
  cat $sheet | sed '1d' | cut -f1 | sort | join -v1 $t - | scripts-send2err
  rm -f $t
  exit 1
endif
rm -f $t

# file names should exist in fastq-or-alignments directory
set fastq_cols = `cat $sheet | head -1 | tr '\t' '\n' | grep -n '^fastq-r[12]$' | cut -d':' -f1`
set fastq_cols = `echo $fastq_cols | tr ' ' ','`
set files = `cat $sheet | sed '1d'  | cut -f$fastq_cols | tr ',\t' '\n' | sort -u | grep -v '^NA$'`
foreach f ($files)
  if (! -e inputs/fastq/$f) then
    scripts-send2err "Error: file \'$f\' not found in inputs/fastq-or-alignments directory."
    exit 1
  endif
end



