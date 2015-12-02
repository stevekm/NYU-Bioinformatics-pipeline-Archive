#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-annotations.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running annotations ============="
resources=1,10G
method=by-object
inpdirs=$(find matrix-filtered/results/ -name '*res_10kb*')                                                                           # run only at 10kb resolution
scripts-master-loop.sh $resources $method ./code/hicseq-annotations.tcsh results/annotations "params/params.*.tcsh" "$inpdirs"


