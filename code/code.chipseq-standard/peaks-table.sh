#!/bin/bash

##
## USAGE: peaks-table.sh OUTPUT-DIR PARAMETER-SCRIPT PEAK-DIRECTORIES
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-bashrc

# process command-line inputs
if (($# != 3)); then
  grep '^##' $0
  exit
fi

# inputs
outdir=$1
params=$2
peakdirs=($3)

# parameters
scripts-send2err "Processing parameter file $params..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- max peak distance = $peakdist"
scripts-send2err "- annotation = $annot_params"

# create merged peaks reference
peaks=( $(echo ${peakdirs[*]} | tr ' ' '\n' | sed 's/$/\/peaks.bed/') )
cat ${peaks[*]} | genomic_regions bed | scripts-sortbed | genomic_regions link -d $peakdist | genomic_regions reg | cut -f2 | rows -number -pref MERGED-PEAK- | sort > $outdir/ref.reg 

# create table
peaks_commas=$(echo ${peaks[*]} | tr ' ' ',')
peaks_labels=( $(echo ${peakdirs[*]} | tr ' ' '\n' | sed 's/.*peaks\.//') )
(echo "MERGED-PEAK-ID ${peaks_labels[*]}" | tr ' ' '\t' ; gtools_threaded matrix -i -nbins 1 --overlap-op hits -format '%.0f' $peaks_commas $outdir/ref.reg | scripts-skipn 1) > $outdir/ref.table.tsv

# annotate merged peaks
(echo "MERGED-PEAK-ID LOCUS PEAK-SIZE REGION GENE-SYMBOL DISTANCE" | tr ' ' '\t'; cat $outdir/ref.reg | genomic_overlaps annotate $annot_params | cut -f1,2,3,4,7 | sed 's/:[^:]\+|/\t/' | sort -u) > $outdir/ref.annotated.tsv

# join tables
scripts-join-tables $outdir/ref.annotated.tsv $outdir/ref.table.tsv > $outdir/peaktable.tsv


