#!/bin/bash
source ./code/code.main/custom-bashrc     # shell settings

# process command-line inputs
if [ $# != 0 ]; then
  echo "##"
  echo "## USAGE: run-domains.sh"
  echo "##"
  exit
fi

# create results directory
scripts-create-path results/

# filter
scripts-send2err "=== Running domains ============="
threads=1,20G
method=by-object
inpdirs=$(find matrix-*/results/ -name '*res_40kb*')                                                                           # run only at 40kb resolution
scripts-master-loop.sh $threads $method ./code/hicseq-domains.tcsh results/domains "params/params.*.tcsh" "$inpdirs"


