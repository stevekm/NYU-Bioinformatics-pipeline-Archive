
module load macs/2.0.10.20131216
set macs_params = "--nomodel --extsize=400"
set use_input = 'true'
set annot_params = "annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 inputs/release/gene-name.bed"

