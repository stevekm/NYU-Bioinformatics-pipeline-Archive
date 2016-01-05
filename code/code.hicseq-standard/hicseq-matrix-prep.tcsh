#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-scale-prep.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
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

# setup
set inpdir = $branch/$object
set filter_branch = ../filter/results/`echo $branch | sed 's/.*\/\(filter\.\)/\1/'`
set n_reads = `cat $filter_branch/$object/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
set res = `echo $bin_size/1000 | bc`
set features = $genome_dir/features-hicnorm/$enzyme.w=${res}kb.features.txt 

# run matrix scaling/normalizing operation
set jid =
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vwE "$chrom_excluded"`)
  scripts-send2err "Processing $mat [enzyme = $enzyme\; features = $features\; n_reads = $n_reads]..."
  
  # estimate required memory
  set mem = `./code/calc-matrix-memory.tcsh $inpdir/$mat 1 5`
  scripts-send2err "requested memory = $mem"

  # run matrix prep
  set jpref = $outdir/__jdata/job.`echo $mat | sed 's/\.[^.]\+$//'`
  scripts-create-path $jpref
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem Rscript ./code/hic-matrix.r normalize -v -o $outdir/$mat --n-reads=$n_reads --ignored-loci=$inpdir/ignored_loci.txt --features=$features $prep_params $inpdir/$mat`)
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


