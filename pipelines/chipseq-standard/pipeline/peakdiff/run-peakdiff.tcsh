#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-peakdiff.tcsh [--dry-run]
##

# ~~~ Entries for auto-report ~~~ #
#DESCRIPTION:
#TITLE: PeakDiff
#FIGURE:
#PARAMS:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set opt = "$1"

# setup
set op = peakdiff
set inpdirs = "inpdirs/matrices"
set filter = "/matrices.[^/]*.nbins_1/"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 1
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v -F "$filter" --exclude-obj="input" "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "group" 2

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



