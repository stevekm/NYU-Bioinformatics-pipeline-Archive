#!/bin/tcsh

##
## USAGE: chipseq-peaks.tcsh OUTPUT-DIR PARAMETER-SCRIPT ALIGNMENT-BRANCH SAMPLES
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

# process command-line inputs
if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0
  exit
endif

set out = $1
set param_script = $2
set branch = $3
set samples = ($4)

# create record of input branch(es) [TODO: need to figure out a way to do a relative symlink]
echo $branch >! $out/obj.branch

# determine input files
set sheet = inputs/sample-sheet.tsv
set treatment_aln = `echo $samples | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`
set control_aln = `echo $samples | tr ' ' '\n' | sed 's/^/^/' | sed 's/$/\t/' | grep -f - $sheet | cut -f3 | grep -vi '^n/a$' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`

# set parameters
source $param_script
scripts-send2err "-- Parameters: "
scripts-send2err "- macs = $macs_params"
scripts-send2err "- annotation = $annot_params"

# run macs2
if ($#control_aln > 0) set control_aln = "-c $control_aln"
set genome = `readlink -f inputs/release/.. | sed 's/.*\///' | cut -d'_' -f1`
set gsize = `grep "^$genome	" params/gsize.info.txt | cut -f2`
if ($gsize == "") then
  scripts-send2err "Error: genome effective size not found in params/gsize.info.txt!"
  exit
endif
scripts-send2err "- genome = $genome"
macs2 callpeak -t $treatment_aln $control_aln --outdir=$out --name=macs $macs_params -g $gsize

# create standardized output files
set ext = narrowPeak
if (`echo $macs_params | grep -cw '\--broad'` == 1) set ext = broadPeak
cat $out/macs_peaks.$ext | cut -f-3,7 >! $out/peak-scores.bed
(echo "PEAK-ID\tCHROMOSOME\tSTART\tEND\tFOLD-CHANGE\tPVALUE(-10log10)\tQVALUE(-10log10)"; cat $out/macs_peaks.$ext | sort -k7,7rg | cut -f-3,7- | rows -number -pref PEAK_) >! $out/peaks.table.tsv
cat $out/peaks.table.tsv | scripts-skipn 1 | cols -t 1 2 3 0 4 | sort -u | sort -k1,1 -k2,2g >! $out/peaks.bed

# annotate peaks
(echo "MERGED-PEAK-ID LOCUS PEAK-SIZE GENE-SYMBOL REGION DISTANCE-START DISTANCE-END" | tr ' ' '\t'; cat $out/peaks.bed | genomic_overlaps $annot_params | cut -f-7 | sort -u) >! $out/peaks.annotated.tsv
scripts-join-tables $out/peaks.table.tsv $out/peaks.annotated.tsv >! $out/peaks.tsv

rm -f $out/peaks.table.tsv $out/peaks.annotated.tsv

