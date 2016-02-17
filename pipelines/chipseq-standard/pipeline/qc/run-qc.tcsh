#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-qc.tcsh [--dry-run]
##

# ~~~ Entries for auto-report ~~~ #
#TITLE: Quality Control
#DESCRIPTION: Check the samples for quality control analysis, produce fingerprint plots.
#FIGURE: chip-fingerprint.pdf
#PARAMS: params.default.tcsh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set opt = "$1"

# setup
set op = fingerprint
set inpdirs = "inpdirs/*"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 4
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "group" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



