#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # custom shell environment

##
## USAGE: read-sample-sheet.tcsh SAMPLE-SHEET-TSV SAMPLE-NAME FEATURE-NAME
##

if ($#argv != 3) then
  grep '^##' $0
  exit 1
endif

set sheet = $1
set sample = $2
set feature = $3

set col = `head -1 $sheet | tr '\t' '\n' | grep -n "^$feature"'$' | cut -d':' -f1`

cat $sheet | awk -F '\t' -v s=$sample '$1==s' | cut -f$col





