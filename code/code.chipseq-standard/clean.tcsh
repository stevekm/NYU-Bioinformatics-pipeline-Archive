#!/bin/tcsh

##
## USAGE: clean
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

rm -rf results
rm -f _tmp*

