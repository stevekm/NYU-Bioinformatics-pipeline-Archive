#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-filtered.tcsh OUTPUT-DIR PARAM-SCRIPT FILTER-BRANCH OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir enzyme"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "params = $matrix_params"

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# setup
./code/create-ignored-loci.tcsh $genome_dir $bin_size >! $outdir/ignored_loci.txt

# generate matrix
set filtered_reads = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/filtered.reg.gz"}'`
scripts-smartcat $filtered_reads | gtools-hic matrix -v -p $outdir/matrix. $matrix_params

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."

