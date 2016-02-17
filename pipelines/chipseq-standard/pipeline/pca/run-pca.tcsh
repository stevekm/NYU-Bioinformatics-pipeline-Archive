#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-pca.tcsh [--dry-run]
##

# ~~~ Entries for auto-report ~~~ #
#DESCRIPTION: This step of the pipeline performs PCA.
#TITLE: PCA
#FIGURE: report_qnorm_page3.pdf
#PARAMS: params.standard.tcsh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set opt = "$1"

# setup
set op = pca
set inpdirs = "inpdirs/matrices"
set filter = "/matrices.[^/]*.nbins_1/"
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 1
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/chipseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v --filter-branch="$filter" "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "*" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



