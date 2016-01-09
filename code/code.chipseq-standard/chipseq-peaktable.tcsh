#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-peaktable.tcsh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH [OBJECTS]
##

# process command-line inputs
if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

# inputs
set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if samples is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- max peak distance = $peakdist"
scripts-send2err "- annotation = $annot_params"

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# create merged peaks reference
set peaks = $branch/*/peaks.bed
cat $peaks | gtools-regions bed | scripts-sortbed | gtools-regions link -d $peakdist | gtools-regions reg | cut -f2 | tools-rows -number -pref MERGED-PEAK- | sort >! $outdir/ref.reg 

# create table
set peaks_commas = `echo $peaks | tr ' ' ','`
set peaks_labels = `echo $peaks | tr ' ' '\n' | sed 's/\/peaks.bed$//' | sed 's/.*\///' | tr '\n' ' '`
(echo "MERGED-PEAK-ID $peaks_labels" | tr ' ' '\t' ; gtools-threaded matrix -i -nbins 1 --overlap-op hits -format '%.0f' $peaks_commas $outdir/ref.reg | scripts-skipn 1) >! $outdir/ref.table.tsv

# annotate merged peaks
(echo "MERGED-PEAK-ID LOCUS PEAK-SIZE GENE-SYMBOL REGION DISTANCE-START DISTANCE-END" | tr ' ' '\t'; cat $outdir/ref.reg | gtools-overlaps $annot_params | cut -f-7 | sort -u) >! $outdir/ref.annotated.tsv

# join tables
scripts-join-tables $outdir/ref.annotated.tsv $outdir/ref.table.tsv >! $outdir/peaktable.tsv

# cleanup
rm -f $outdir/ref.annotated.tsv $outdir/ref.table.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


