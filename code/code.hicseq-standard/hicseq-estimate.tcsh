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
cat $genome_dir/genome.bed | genomic_regions win -s $bin_size -d $bin_size | genomic_overlaps subset -i $genome_dir/centrotelo.bed | genomic_regions shiftp -5p 1 -3p 0 | sed 's/\t/:/' | sed 's/\t/-/' >! $outdir/ignored_loci.txt

# run estimation
set inpdir = $branch/$sample
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM'`)
  send2err "Processing input matrix $mat..."
  foreach g ($gamma)
    set outmat = `echo $mat | sed 's/.tsv$/'".gamma=$g.RData/"`
    set jpref = $outdir/job.$outmat
    set jid = `scripts-qsub-run $jpref 1 40 ./code/hic_matrix.r estimate -v -o $outdir/$outmat --ignored-loci=$outdir/ignored_loci.txt --gamma=$g $hic_params $inpdir/$mat`
    scripts-qsub-wait $jid
    set t = `scripts-create-temp $outdir`
    tail $jpref.out >! $t
    cat $t >! $jpref.out
    rm -f $t
  end
end

# save variables
set >! $outdir/job.vars.tsv


