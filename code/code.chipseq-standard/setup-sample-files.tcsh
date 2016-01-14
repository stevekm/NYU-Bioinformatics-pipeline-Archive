#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: setup-sample-files.tcsh
##

# process command-line inputs
if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# determine sample names
set inpdir = fastq.raw
set outdir = fastq
set samples = `cd $inpdir; ls -1d *.fastq.gz *.bam | sed 's/\.fastq.gz$//' | sed 's/\.bam$//'`      # TODO: handle R1/R2, multiple lanes etc
foreach sample ($samples)
  mkdir $outdir/$sample
  foreach f ($inpdir/${sample}*.fastq.gz $inpdir/${sample}*.bam) 
    (cd $outdir/$sample; ln -s ../../$f)
  end
end


