#!/bin/tcsh

##
## USAGE: easydiff-filter.tcsh EASYDIFF-DIRECTORY PARAMETERS
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 2) then
  grep '^##' $0
  exit
endif

set d = $1
set params = "$2"

eval $params
if (! $?min_log2fc) set min_log2fc = 1.0
if (! $?min_val) set min_val = 2.0

head -1 $d/diff.score
cat $d/diff.score | scripts-skipn 1 | awk -v x1=$min_log2fc -v x2=$min_val '$2>=x1 && $7>=x2' | sort -k7,7rg
cat $d/diff.score | scripts-skipn 1 | awk -v x1=$min_log2fc -v x2=$min_val '$2<=-x1 && $6>=x2' | sort -k6,6rg


