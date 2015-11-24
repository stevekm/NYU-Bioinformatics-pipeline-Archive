#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-ic.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH SAMPLE
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set sample = $4

# run parameter script
source $params

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch/$sample "genome genome_dir enzyme bin_size"

# create path
scripts-create-path $outdir/

# run iteractive correction
set x = `scripts-create-temp $outdir`
set y = `scripts-create-temp $outdir`
set inpdir = $branch/$sample
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM'`)
  scripts-send2err "Processing $mat..."
  cat $inpdir/$mat | scripts-skipn 1 | cut -f2- >! $x
  python ./code/hicseq-ic.py $x $cutoff $y
  cat $inpdir/$mat | head -1 >! $outdir/$mat
  cat $inpdir/$mat | scripts-skipn 1 | cut -f1 | paste - $y | matrix format -n 3 | sed 's/ *$//' | tr -s ' ' '\t' >> $outdir/$mat
end

# cleanup
rm -f $x $y

# save variables
set >! $outdir/job.vars.tsv


