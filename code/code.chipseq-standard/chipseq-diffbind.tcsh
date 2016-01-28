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

# setup sample sheet for diffbind
set diffbind_sample_sheet = $outdir/diffbind-sample-sheet.csv
echo "header-info..." >! $diffbind_sample_sheet
foreach obj ($objects)
  set treatment_bam = `cat $branch/$obj/job.vars.tsv | grep "^treatment_aln	" | cut -f2`
  set control_bam = `cat $branch/$obj/job.vars.tsv | grep "^control_aln	" | cut -f2 | sed 's/^.* //'`
  set macs_peaks_xls = $branch/$obj/macs_peaks.xls
  set obj_group = `./code/read-sample-sheet.tcsh $sheet $obj $group_var`
  echo "$obj,$treatment_bam,$control_bam,$macs_peaks_xls,$obj_group" >> $diffbind_sample_sheet
end

# run diffbind
# TODO: remember to use diffbind_param variable (which is set in the params.standard.tcsh script), instead of hard-coding the params. 


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




