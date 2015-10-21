#!/bin/tcsh

##
## USAGE: validate-sample-sheet SAMPLE-SHEET-TSV
##

if ($#argv != 1) then
  grep '^##' $0
  exit
endif

set sheet = $1

# set path
set path = (./code/code $path)

if (`cat $sheet | grep -v '^#' | grep ' ' | wc -l` > 0) then
  scripts-send2err 'Error: spaces are not allowed in sample sheet.'
  exit
endif

if (`cat $sheet | grep -v '^#' | cut -f2 | sort | uniq -d | wc -l` > 0) then
  scripts-send2err 'Error: duplicate sample names are not allowed in sample sheet.'
  exit
endif




