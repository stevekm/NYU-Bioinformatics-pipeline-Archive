#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-heatmaps.tcsh [--dry-run]
##

# ~~~ Entries for auto-report ~~~ #
#TITLE: Heatmaps
#XFIGUREX: clustering.png #Keep this disabled until convert code is added
#PARAMS: params.mean.tcsh params.standard.tcsh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set opt = "$1"

# setup
set op = heatmaps
set inpdirs = "inpdirs/matrices"
set filter = "/matrices.[^/]*.nbins_[0-9][0-9]"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 1
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v --filter-branch="$filter" "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "*" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



