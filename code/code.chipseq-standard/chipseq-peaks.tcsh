#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # shell settings (must be included in all scripts)

##
## USAGE: chipseq-peaks.tcsh OUTPUT-DIR PARAMETER-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- macs = $macs_params"
scripts-send2err "- annotation = $annot_params"

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine input files
set treatment_aln = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`
set control_aln = 
if ($use_input != 'false') then
  set control_samples = `./code/read-sample-sheet.tcsh $sheet "$objects" control | sort -u | grep -v '^NA$'`
  if ("$control_samples" != "") then
    set control_aln = `echo $control_samples | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`
    set control_aln = "-c $control_aln"
  endif
endif

# run macs2
set g = `echo $genome | sed 's/[0-9]\+$//'`
set gsize = `grep "^$g	" params/gsize.info.txt | cut -f2`
if ($gsize == "") then
  scripts-send2err "Error: genome effective size not found in params/gsize.info.txt!"
  exit 1
endif
scripts-send2err "- genome = $genome"
echo macs2 callpeak -t $treatment_aln $control_aln --outdir=$outdir --name=macs $macs_params -g $gsize | scripts-send2err           # log info
macs2 callpeak -t $treatment_aln $control_aln --outdir=$outdir --name=macs $macs_params -g $gsize

# create standardized output files
set ext = narrowPeak
if (`echo $macs_params | grep -cw '\--broad'` == 1) set ext = broadPeak
cat $outdir/macs_peaks.$ext | cut -f-3,7 >! $outdir/peak-scores.bed
(echo "PEAK-ID\tCHROMOSOME\tSTART\tEND\tFOLD-CHANGE\tPVALUE(-10log10)\tQVALUE(-10log10)"; cat $outdir/macs_peaks.$ext | sort -k7,7rg | cut -f-3,7- | tools-rows -number -pref PEAK_) >! $outdir/peaks.table.tsv
cat $outdir/peaks.table.tsv | scripts-skipn 1 | tools-cols -t 1 2 3 0 4 | sort -u | sort -k1,1 -k2,2g >! $outdir/peaks.bed

# annotate peaks
(echo "MERGED-PEAK-ID LOCUS PEAK-SIZE GENE-SYMBOL REGION DISTANCE-START DISTANCE-END" | tr ' ' '\t'; cat $outdir/peaks.bed | gtools-overlaps $annot_params | cut -f-7 | sort -u) >! $outdir/peaks.annotated.tsv
scripts-join-tables $outdir/peaks.table.tsv $outdir/peaks.annotated.tsv >! $outdir/peaks.tsv

# cleanup
rm -f $outdir/peaks.table.tsv $outdir/peaks.annotated.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



