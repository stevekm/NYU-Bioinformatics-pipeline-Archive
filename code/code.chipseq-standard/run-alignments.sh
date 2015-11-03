#!/bin/bash

##
## USAGE: run-alignments.sh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]
then
  grep '^##' $0
  exit
fi

# parameters
params=params/params.bowtie2.tcsh         # parameter script needs to be written in same language as the wrapper script!

# create results directory
scripts-create-path results/

# read from sample sheet(s)
sheet=inputs/sample-sheet.tsv
./code/validate-sample-sheet.tcsh $sheet
Q=( $(cat $sheet | grep -v '^#' | cut -f4) )           # file names (fastq or bam)
S=( $(cat $sheet | grep -v '^#' | cut -f1) )           # sample names

# align
scripts-send2err "=== Aligning reads ============="
threads=8
jid=()
k=0
while (( $k < ${#Q[@]} )); do
  q=$(echo ${Q[$k]} | tr ',' '\n' | sed 's/^/inputs\/fastq-or-alignments\//' | tr '\n' ',' | sed 's/,$//')      # this allows comma-separated fastq files
  qname=${S[$k]}
  outdir=results/$qname
  jid+=( $(scripts-qsub-wrapper $threads ./code/chipseq-align.tcsh $outdir $params $q) )
  let k=k+1
done

# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



