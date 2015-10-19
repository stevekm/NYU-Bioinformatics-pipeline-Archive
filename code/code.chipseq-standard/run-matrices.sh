#!/bin/bash

##
## USAGE: run-matrices.sh
##

if [ $# != 0 ]; then
  grep '^##' $0
  exit
fi

# set path
PATH=./code/code:$PATH

# generate matrices
scripts-send2err "=== Generating matrices ============="
scripts-create-path results/
set jid = ()
foreach params (params/params.*)
  # parameters
  scripts-send2err "Processing parameter file $params..."
  source $params
  scripts-send2err "- window = $win"
  scripts-send2err "- flank = $flank"
  scripts-send2err "- nbins = $nbins"
 
  # loop over all peak parameter settings
  set PEAK_PARAMS = `find peaks/results/ -name 'peaks.bed' | sed 's/\/peaks\.[^/]\+\/peaks.bed//' | sort -u`
  foreach peak_params ($PEAK_PARAMS)
    # create output directory
    set out = results/matrices.`echo $params | sed 's/.*\///' | sed 's/^params\.//'`/`echo $peak_params | sed 's/peaks\/results\///'`
    set regions = $peak_params/peaks.*/peaks.bed
    scripts-send2err "-- out = $out"
    
    # create ref.bed
    set ref = $out/ref.bed
    if (! -e $ref) then
      scripts-create-path $ref
      if ($win > 0) then
        set p = `scripts-create-temp`
        cat $regions | genomic_regions bed >! $p
        cat inputs/release/genome.bed | genomic_regions win -s $win -d `echo $win/4 | bc` | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | genomic_overlaps subset -i $p >! $ref
        rm -f $p
      else
        cat $regions | genomic_regions bed | scripts-sortbed | genomic_regions link | genomic_regions center | genomic_regions pos -op 1 | genomic_regions shiftp -5p -$flank -3p +$flank | genomic_regions shiftp -5p 0 -3p -1 | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >! $ref
      endif
    endif
  
    # create matrices
    scripts-create-path $out/logs/
    set nthreads = 2
    foreach aln (alignments/results/align.*)
      set sample = `echo $aln | sed 's/.*\/align\.//'`
      set bam = $aln/alignments.bam
      set mat = $out/matrix.$sample.tsv
      if (! -e $mat) then
        scripts-send2err "Processing $sample..." 
        scripts-create-path $mat
        set jid = ($jid `scripts-qsub-run $out/logs/job.matrix.$sample $nthreads ./code/create-matrix $mat $nthreads $bam $ref $nbins`)
      else 
        scripts-send2err "Warning: $mat already exists, skipping..."
      endif
    end
  end
end

scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "$jid"
scripts-send2err "Done."



