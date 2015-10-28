#!/bin/bash

##
## USAGE: run-align.sh
##

# shell settings
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]
then
  grep '^##' $0
  exit
fi

# create results directory
scripts-create-path results/

# read from sample sheet(s)
sheet=inputs/sample-sheet.tsv
#./code/validate-sample-sheet.tcsh $sheet
S=( $(cat $sheet | grep -v '^#' | cut -f1) )            # sample names
Q1=( $(cat $sheet | grep -v '^#' | cut -f2) )           # fastq R1 files
Q2=( $(cat $sheet | grep -v '^#' | cut -f3) )           # fastq R2 files

# align
scripts-send2err "=== Aligning read pairs ============="
threads=8
jid=()
# loop over all parameter files
for params in $(ls -1 params/params.*.tcsh); do
  pname=$(echo $params | sed 's/.*\///' | cut -d'.' -f2)

  # loop over all samples
  k=0
  while (( $k < ${#S[@]} )); do
    q1=$(echo ${Q1[$k]} | tr ',' '\n' | sed 's/^/inputs\/fastq\//')      # allows comma-separated R1 fastq files
    q2=$(echo ${Q2[$k]} | tr ',' '\n' | sed 's/^/inputs\/fastq\//')      # allows comma-separated R2 fastq files
    qname=${S[$k]}
    outdir=results/align.$pname/$qname
    jid+=( $(scripts-qsub-wrapper $threads ./code/hicseq-align.tcsh $outdir $params "$q1" "$q2") )
    let k=k+1
  done
done
       
# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



