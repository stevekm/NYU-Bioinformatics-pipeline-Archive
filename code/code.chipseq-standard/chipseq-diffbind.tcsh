#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-diffbind.tcsh OUTPUT-DIR PARAM-SCRIPT PEAKS-BRANCH [OBJECTS]
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

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

echo "sheet: $sheet"
echo "genome: $genome"
echo "branch: $branch"
#echo "group_var: $group_var"
echo "diffbind_factor: $diffbind_factor"
# optional blocking factor (for paired analysis)
echo "diffbind_blocking_factor: $diffbind_blocking_factor"

# setup sample sheet for diffbind
set diffbind_sample_sheet = $outdir/diffbind-sample-sheet.csv
echo "SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,bamControl,Peaks,PeakCaller" >! $diffbind_sample_sheet
foreach obj ($objects)
  set treatment_bam = `cat $branch/$obj/job.vars.tsv | grep "^treatment_aln	" | cut -f2`
  set control_bam = `cat $branch/$obj/job.vars.tsv | grep "^control_aln	" | cut -f2 | sed 's/^.* //'`
  set macs_peaks_xls = $branch/$obj/macs_peaks.xls
  set db_factor = `./code/read-sample-sheet.tcsh $sheet $obj $diffbind_factor`
  if ("$diffbind_blocking_factor" == "") then
    set db_block_factor = ""
  else
    set db_block_factor = `./code/read-sample-sheet.tcsh $sheet $obj $diffbind_blocking_factor`
  endif
  echo "$obj,?,?,$db_factor,?,$db_block_factor,$treatment_bam,$control_bam,$macs_peaks_xls,macs" >> $diffbind_sample_sheet
end

Rscript --vanilla ./code/chipseq-diffbind.R $outdir $diffbind_sample_sheet $genome $diffbind_blocking_factor


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




