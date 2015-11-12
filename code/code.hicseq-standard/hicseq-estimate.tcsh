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

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# setup
set release = inputs/release
cat $release/genome.bed | genomic_regions win -s ${res}000 -d ${res}000 | genomic_overlaps subset -i $release/../DNA/centrotelo.bed | genomic_regions shiftp -5p 1 -3p 0 | sed 's/\t/:/' | sed 's/\t/-/' >! $outdir/ignored_loci.txt

# run estimation
set inpdir = $branch/$sample
foreach mat (`cd $inpdir; ls -1 matrix.*.tsv | grep -vw 'chrM'`)
  send2err "Processing input matrix $mat..."
  foreach g ($gamma)
    set outmat = `echo $mat | sed 's/.tsv$/'".gamma=$g.RData/"`
    ./code/hic_matrix.r estimate -v -o $outdir/$outmat --ignored-loci=$outdir/ignored_loci.txt --gamma=$g $hic_params $inpdir/$mat          # TODO: qsub this, wait for all jobs to finish!
  end
end



