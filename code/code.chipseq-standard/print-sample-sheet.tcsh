#!/bin/tcsh

##
## USAGE: print-sample-sheet SAMPLE-SHEET-TSV
##

if ($#argv != 1) then
  grep '^##' $0
  exit
endif

set sheet = $1

cat $sheet | awk '{ printf "%-50s %-40s %-40s %-20s\n", $1, $2, $3, $4 }'


