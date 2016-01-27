#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: chipseq-fingerprint.tcsh OUTPUT-DIR PARAMETER-SCRIPT ALIGNMENT-BRANCH OBJECTS
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0 | scripts-send2err
  exit
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

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine input files
set control_objects = `./code/read-sample-sheet.tcsh $sheet "$objects" control | sort -u | grep -v '^NA$'`
set labels = `echo $objects $control_objects | tr ' ' '\n' | sort -u`
set bam = `echo $labels | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`

# run deeptools wrapper
echo ./code/chipseq-fingerprint-deeptools.tcsh $outdir $params "'$bam'" "'$labels'" | scripts-send2err 
./code/chipseq-fingerprint-deeptools.tcsh $outdir $params "$bam" "$labels"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



