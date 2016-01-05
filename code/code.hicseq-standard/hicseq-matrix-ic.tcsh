#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-ic.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT(S)
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

# run iteractive correction
set inpdir = $branch/$object
set jid = 
set matrices = `cd $inpdir; ls -1 matrix.*.tsv | grep -vwE "$chrom_excluded"`
foreach mat ($matrices)
  scripts-send2err "Processing $mat..."
  
  # estimate required memory
  set mem = `./code/calc-matrix-memory.tcsh $inpdir/$mat 1.5 5`
  scripts-send2err "requested memory = $mem"

  # run IC
  set jpref = $outdir/__jdata/job.`echo $mat | sed 's/\.[^.]\+$//'`
  scripts-create-path $jpref
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem ./code/hicseq-matrix-ic-mirny.tcsh $outdir $params $branch $object $mat`)
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


