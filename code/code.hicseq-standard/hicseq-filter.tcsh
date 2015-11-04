#!/bin/tcsh

##
## USAGE: hic-filter.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-BRANCH SAMPLE
##

# shell settings
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set aln_branch = $3
set sample = $4

# run parameter script
source $params

# create path
scripts-create-path $out/

# determine enzyme name
if (`echo $sample | grep -ic '\-NcoI'` == 1) then
  set enzyme = NcoI
else if (`echo $sample | grep -ic '\-HindIII'` == 1) then
  set enzyme = HindIII
else
  scripts-send2err "Error: enzyme information not embedded in file name, aborting..."
  exit
endif

# filter
samtools view $aln_branch/$sample/alignments.bam | gtools_hic filter -v -E $release/../DNA/$enzyme.fragments.bed $filter_params | gzip >! $out/filtered.reg.gz






