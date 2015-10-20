#!/bin/tcsh

##
## USAGE: run.alignments
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# set path
set path = (./code/code $path)

# parameters
set params = params/params.bowtie2.tcsh

# create results directory
scripts-create-path results/

# read from sample sheet(s)
set sheet = inputs/sample-sheet.tsv
if (`cat $sheet | grep -v '^#' | grep ' ' | wc -l` > 0) then
  scripts-send2err 'Error: spaces are not allowed in sample sheet.'
  exit
endif
set Q = `cat $sheet | grep -v '^#' | cut -f1`           # file names (fastq or bam)
set S = `cat $sheet | grep -v '^#' | cut -f2`           # sample names
if (`echo $S | tr ' ' '\n' | sort | uniq -d | wc -l` > 0) then
  scripts-send2err 'Error: duplicate sample names are not allowed.'
  exit
endif

# align
scripts-send2err "=== Aligning reads ============="
set threads = 8
set jid = ()
set k = 1
while ($k <= $#Q)
  set q = `echo $Q[$k] | tr ',' '\n' | sed 's/^/inputs\/fastq-or-alignments\//' | tr '\n' ',' | sed 's/,$//'`   # this allows comma-separated fastq files
  set qname = $S[$k]
  set outdir = results/align.$qname
  set jid = ($jid `scripts-qsub-wrapper $threads ./code/chipseq-align $outdir $params $q`)
  @ k ++
end

# wait for all processes
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "$jid"
scripts-send2err "Done."



