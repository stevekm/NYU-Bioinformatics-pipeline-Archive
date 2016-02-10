#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # custom shell environment

##
## USAGE: read-sample-sheet.tcsh SAMPLE-SHEET-TSV SAMPLE-NAME(S) FEATURE-NAME [INCLUDE-SAMPLE-NAME=no]
##

if ($#argv < 3) then
  grep '^##' $0
  exit 1
endif

set sheet = $1
set samples = "$2"
set feature = $3
set include = $4

set col = `head -1 $sheet | tr '\t' '\n' | grep -n "^$feature"'$' | cut -d':' -f1`
if ("$include" == 'yes') then
  set col = "1,$col"
endif

if ("$samples" == "*") then
  cat $sheet | sed '1d' | cut -f$col
else
  set t = `scripts-create-temp`
  echo $samples | tr ' ' '\n' | sort -u >! $t
  cat $sheet | sed '1d' | sort | join -t '	' $t - | cut -f$col
  rm -f $t
endif






