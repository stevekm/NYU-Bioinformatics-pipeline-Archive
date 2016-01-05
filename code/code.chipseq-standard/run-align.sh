#!/bin/bash
source ./code/code.main/custom-bashrc    # shell settings (must be included in all scripts)

##
## USAGE: run-align.sh
##

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# setup
op=align
params=params/params.bowtie2.tcsh         # parameter script needs to be written in same language as the wrapper script!
scripts-send2err "=== Operation = $op ============="

# create results directory
scripts-create-path results/

# read from sample sheet(s)
sheet=inputs/sample-sheet.tsv
./code/validate-sample-sheet.tcsh $sheet
S=( $(./code/read-sample-sheet.tcsh $sheet '*' sample) )           # sample names

# run
threads=8
jid=()
k=0
for sample in ${S[*]}; do
  outdir=results/$sample
  jid+=( $(scripts-qsub-wrapper $threads ./code/chipseq-$op.tcsh $outdir $params inputs/fastq $sample) )
done

# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"
scripts-send2err "Done."



