#!/bin/bash

##
## USAGE: chipseq-peaktable.sh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH
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
peaks_branch=$3

# create output path
scripts-create-path $outdir/

# create record of input branch(es) [TODO: need to figure out a way to do a relative symlink]
echo $peaks_branch > $outdir/obj.branch

# parameters
scripts-send2err "Processing parameter file $params..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- max peak distance = $peakdist"
scripts-send2err "- annotation = $annot_params"

# create merged peaks reference
peaks=$peaks_branch/*/peaks.bed
cat $peaks | genomic_regions bed | scripts-sortbed | genomic_regions link -d $peakdist | genomic_regions reg | cut -f2 | rows -number -pref MERGED-PEAK- | sort > $outdir/ref.reg 

# create table
peaks_commas=$(echo $peaks | tr ' ' ',')
peaks_labels=$(echo $peaks | tr ' ' '\n' | sed 's/\/peaks.bed$//' | sed 's/.*\///' | tr '\n' ' ')
(echo "MERGED-PEAK-ID $peaks_labels" | tr ' ' '\t' ; gtools_threaded matrix -i -nbins 1 --overlap-op hits -format '%.0f' $peaks_commas $outdir/ref.reg | scripts-skipn 1) > $outdir/ref.table.tsv

# annotate merged peaks
(echo "MERGED-PEAK-ID LOCUS PEAK-SIZE GENE-SYMBOL REGION DISTANCE-START DISTANCE-END" | tr ' ' '\t'; cat $outdir/ref.reg | genomic_overlaps $annot_params | cut -f-7 | sort -u) > $outdir/ref.annotated.tsv

# join tables
scripts-join-tables $outdir/ref.annotated.tsv $outdir/ref.table.tsv > $outdir/peaktable.tsv

# cleanup
rm -f $outdir/ref.annotated.tsv $outdir/ref.table.tsv


