#!/bin/tcsh
source ./code/code.main/custom-tcshrc      # customize shell environment

##
## USAGE: run-domains.tcsh [--dry-run]
##

#% This step identifies topologically-associated domains (TADs) using different methods.   
#% TABLES: 
#% FIGURES:

# process command-line inputs
if ($#argv > 1) then
  grep '^##' $0 | scripts-send2err
  exit
endif

set opt = "$1"

# setup
set op = domains
set inpdirs = "inpdirs/*"
set filter = "*.res_40kb"                  # work only with 40kb resolution
set results = results
scripts-create-path $results/
scripts-send2err "=== Operation = $op ============="
set resources = 1,20G
set cmd = "./code/code.main/scripts-qsub-wrapper $resources ./code/hicseq-$op.tcsh"

# generate run script
Rscript ./code/code.main/pipeline-master-explorer.r -v -F "$filter" --exclude-branch='/matrix-estimated.[^/]+prep_distlog2' "$cmd" $results/$op "params/params.*.tcsh" "$inpdirs" "" "sample" 1

# run and wait until done!
if ("$opt" != "--dry-run") scripts-submit-jobs ./$results/.db/run



