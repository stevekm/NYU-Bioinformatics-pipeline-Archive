#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-estimate.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)          # TODO: allow multiple samples, so that this operation can be run by-group

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# create table of interactions for each chromosome
set inpdir = $branch/$object
set filter_branch = ../filter/results/`echo $branch | sed 's/.*\/\(filter\.\)/\1/'`
set n_reads = `cat $filter_branch/$object/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
set jid = 
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`)
  scripts-send2err "Processing input matrix $inpdir/$mat..."
  set outmatdir = `echo $mat | sed 's/\.tsv$//' | sed 's/\.RData$//'`
  set mem = `./code/calc-matrix-memory.tcsh $inpdir/$mat 2 5`
  scripts-send2err "requested memory = $mem"
  set jdata = $outdir/__jdata
  scripts-create-path $jdata/
  set jpref = $jdata/job.$outmatdir
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem Rscript ./code/hic-matrix.r loops -v -o $outdir/$outmatdir -L $object --n-reads=$n_reads $loop_params $inpdir/$mat`)
end
scripts-qsub-wait "$jid"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


