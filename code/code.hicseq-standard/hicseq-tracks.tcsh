#!/bin/tcsh

source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-tracks.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

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
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Create tracks
set sample_paths = `echo "$branch\t$objects" | tools-key-expand | tr '\t' '/'`
set sample_reads = `echo $sample_paths | tr ' ' '\n' | sed 's/$/\/filtered.reg.gz/'`
./code/hicseq-tracks-washu.tcsh $outdir "$sample_reads" $genome_dir/genome.bed $bin_size

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."

