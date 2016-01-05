#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-interactions.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Identifying interactions ============="
resources=1
op=interactions
method=by-object
inpdirs=$(find matrix-filtered/results/ -name '*res_10kb.maxd*')                                                                           # run only at 10kb resolution
scripts-master-loop.sh $resources $method ./code/hicseq-$op.tcsh results/$op "params/params.*.tcsh" "$inpdirs"


