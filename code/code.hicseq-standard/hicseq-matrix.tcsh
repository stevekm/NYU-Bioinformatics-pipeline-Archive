#!/bin/tcsh

##
## USAGE: hicseq-matrix.tcsh OUTPUT-DIR PARAM-SCRIPT FILTERED-DIR
##

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set filtered = ($3)

# set path
set path = (./code/code $path)

# run parameter script
source $params

# create path
scripts-create-path $out

# generate matrix
set filtered_reads = `echo $filtered | tr ' ' '\n'  | sed 's/$/\/filtered.reg.gz/'`
scripts-smartcat $filtered_reads | gtools_hic matrix -v -p $out/matrix. $matrix_params



