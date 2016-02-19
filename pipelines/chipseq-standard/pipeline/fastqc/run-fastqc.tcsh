#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-fastqc.tcsh [--dry-run]
##

# ~~~ Entries for auto-report ~~~ #
#TITLE: FastQC
#DESCRIPTION: Run FastQC on the input fastq files.
#FIGURE:
#SHEET:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit
endif

set opt = "$1"

# setup
set op = fastqc
set inpdirs = "inputs"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 8 # use a range since we might be calling this on many files at once
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# symlink fastq to results (required for first step of the pipeline)
if (! -e inputs/$results) then
  (cd inputs; ln -s fastq $results)
endif

# ~~ develpoment
echo "> " "$0 OUTPUT HERE"
echo "> " "$cmd IS cmd"
echo "> " "$results/$op IS results/op"
echo "> " "params/params.*.tcsh" "IS PARAMS"
echo "> " "$inpdirs IS inpdirs"
# echo  " ("$cmd") $results/$op "params/params.*.tcsh" ("$inpdirs") "" ("sample") 1 "
# ~~~

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "sample" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run
