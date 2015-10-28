#!/bin/bash

##
## USAGE: run-peaktable.sh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# generate peak tables
func=peaktable
scripts-send2err "=== Generating peak tables ============="
scripts-create-path results/
sheet=inputs/sample-sheet.tsv
threads=1
jid=()
# look over all peaktable params
for params in $(ls -1 params/params.*.sh); do
  params_name=$(echo $params | sed 's/.*\///' | sed 's/^params\.//' | sed 's/\.sh$//')
  
  # loop over all peak params
  for peak_params in $(find peaks/results/ -name 'peaks.bed' | sed 's/\/peaks\.[^/]\+\/peaks.bed//' | sort -u); do
    peak_params_name=$(echo $peak_params | sed 's/peaks\/results\///')

    # generate output directory name
    out=results/$func.$params_name/$peak_params_name
    scripts-send2err "-- out = $out"
    
    # run
    peakdirs=( $(ls -1d $peak_params/peaks.* | grep -vi input) )                             # excluding inputs
    jid+=( $(scripts-qsub-wrapper $threads ./code/peaks-table.sh $out $params "${peakdirs[*]}") )

  done
done

# wait
scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "${jid[*]}"

# done
scripts-send2err "Done."



