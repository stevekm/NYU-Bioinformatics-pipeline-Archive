#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-estimate.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH SAMPLE
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set sample = $4

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch/$sample "genome genome_dir enzyme bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# setup
./code/create-ignored-loci.tcsh $genome_dir $bin_size >! $outdir/ignored_loci.txt

# run estimation
set inpdir = $branch/$sample
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM'`)
  scripts-send2err "Processing input matrix $mat..."
  set n_rows = `cat $inpdir/$mat | wc -l`
  set n_cols = `head -1 $inpdir/$mat | tr '\t' '\n' | wc -l`
  set mem = `echo "500*$n_rows*$n_cols/1000000000+1" | bc`
  scripts-send2err "requested memory = $mem"
  foreach g ($gamma)
    set outmat = `echo $mat | sed 's/.tsv$/'".gamma=$g.RData/"`
    set jpref = $outdir/job.$outmat
    set jid = `scripts-qsub-run $jpref 1 $mem Rscript ./code/hic-matrix.r estimate -v -o $outdir/$outmat --ignored-loci=$outdir/ignored_loci.txt --gamma=$g $hic_params $inpdir/$mat`
    scripts-qsub-wait $jid
    set t = `scripts-create-temp $outdir`
    tail $jpref.out >! $t
    cat $t >! $jpref.out
    rm -f $t
  end
end

# save variables
set >! $outdir/job.vars.tsv


