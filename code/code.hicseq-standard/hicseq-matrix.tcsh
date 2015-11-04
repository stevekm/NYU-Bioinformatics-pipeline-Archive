#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix.tcsh OUTPUT-DIR PARAM-SCRIPT FILTER-BRANCH SAMPLES
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set flt_branch = $3
set samples = ($4)

# run parameter script
source $params

# create path
scripts-create-path $out

# generate matrix
set filtered_reads = `echo $samples | tr ' ' '\n' | awk -v d=$flt_branch '{print d"/"$0"/filtered.reg.gz"}'`
scripts-smartcat $filtered_reads | gtools_hic matrix -v -p $out/matrix. $matrix_params



