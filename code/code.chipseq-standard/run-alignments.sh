#!/bin/bash

##
## USAGE: run-alignments.sh
##

if [ $# != 0 ]
then
  grep '^##' $0
  exit
fi

# set path
PATH=./code/code:$PATH

# parameters
params=params/params.bowtie2.tcsh

# create results directory
scripts-create-path results/

# read from sample sheet(s)
sheet=inputs/sample-sheet.tsv
if [ $(cat $sheet | grep -v '^#' | grep ' ' | wc -l) -gt 0 ]; then
  scripts-send2err 'Error: spaces are not allowed in sample sheet.'
  exit
fi
Q=( $(cat $sheet | grep -v '^#' | cut -f1) )           # file names (fastq or bam)
S=( $(cat $sheet | grep -v '^#' | cut -f2) )           # sample names
if [ $(echo $S | tr ' ' '\n' | sort | uniq -d | wc -l) -gt 0 ]; then
  scripts-send2err 'Error: duplicate sample names are not allowed.'
  exit
fi

# align
scripts-send2err "=== Aligning reads ============="
threads=8
jid=()
k=0
while [ $k -lt ${#Q[@]} ]; do
  q=$(echo ${Q[$k]} | tr ',' '\n' | sed 's/^/inputs\/fastq-or-alignments\//' | tr '\n' ',' | sed 's/,$//')      # this allows comma-separated fastq files
  qname=${S[$k]}
  outdir=results/align.$qname
  jid+=( $(scripts-qsub-wrapper $threads ./code/chipseq-align $outdir $params $q) )
  let k=k+1
done
         
# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



