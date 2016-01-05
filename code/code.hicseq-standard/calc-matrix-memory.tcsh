#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: calc-matrix-memory.tcsh MATRIX FACTOR MIN-MEMORY [ZONE-SIZE]
##

# process command-line inputs
if ($#argv < 3) then
  grep '^##' $0
  exit
endif

set mat = $1
set n = $2
set min_mem = $3
set zone_size = $4

set n_rows = `cat $mat | wc -l`
if (($zone_size == "") || ($zone_size == 0)) then
  set n_cols = `head -1 $mat | tr '\t' '\n' | wc -l`
else
  set n_cols = $zone_size
endif
set mem = `echo "100*$n*$n_rows*$n_cols/1000000000" | bc`         # in Gb
if ($mem <= $min_mem) set mem = $min_mem

echo ${mem}G


