#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc
module load macs/2.0.10.20131216

set extsize = `./code/read-sample-sheet.tcsh $sheet "$objects" fragmentation-size`
set extsize = `echo $extsize | tools-vectors m -n 0`
set macs_params = "--broad --nomodel --extsize=$extsize"
set use_input = 'true'
set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 $genome_dir/gene-name.bed"

