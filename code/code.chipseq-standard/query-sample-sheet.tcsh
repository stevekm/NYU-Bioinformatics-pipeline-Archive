#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # custom shell environment

##
## USAGE: query-sample-sheet.tcsh FIELD SAMPLE-NAMES
##

if ($#argv != 2) then
  grep '^##' $0
  exit 1
endif

set field = $1
set samples = ($2)

set sheet = inputs/sample-sheet.tsv

if ($field == 'control') then
  echo $samples | tr ' ' '\n' | sed 's/^/^/' | sed 's/$/\t/' | grep -f - $sheet | cut -f3 | grep -vi '^n/a$'

else
  scripts-send2err "Error: [query-sample-sheet.tsch] field name does not exist."
  exit 1
endif


