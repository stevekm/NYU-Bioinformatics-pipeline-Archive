#!/bin/tcsh

##
## USAGE: create-matrix OUTPUT-DIR N-THREADS ALIGNED-READS REF-REGION NBINS [OVERLAP-OP={hits,value}]
##

if ($#argv < 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set n_proc = $2
set reads = ($3)
set ref = $4
set nbins = $5
set op = $6
if ($op == '') set op = hits

# set path
set path = (./code/code $path)

set reads = `echo $reads | tr ' ' ','`
if (`genomic_regions reg $ref | cut -f1 | sort | uniq -d | wc -l` > 0) then
  scripts-send2err "Error: unique labels required in the $ref file."
  exit
else
  gtools_threaded matrix -v -i -p $n_proc --overlap-op $op -nbins $nbins -rpkm $reads $ref >! $outdir/matrix.tsv
endif




