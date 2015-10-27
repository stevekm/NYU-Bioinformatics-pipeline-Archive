#!/bin/tcsh

##
## USAGE: hic-filter.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-DIR
##

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set aln = $3

# set path
set path = (./code/code $path)

# run parameter script
source $params

# create path
scripts-create-path $out

# determine enzyme name
if (`echo $aln | grep -ic '\-NcoI'` == 1) then
  set fragments = $release/../DNA/NcoI.fragments.bed
else if (`echo $aln | grep -ic '\-HindIII'` == 1) then
  set fragments = $release/../DNA/HindIII.fragments.bed
else
  scripts-send2err "Error: enzyme information not embedded in file name, aborting..."
  exit
endif

# filter
samtools view $aln/alignments.bam | gtools_hic filter -v -E $fragments $filter_params | gzip >! $out/filtered.reg.gz






