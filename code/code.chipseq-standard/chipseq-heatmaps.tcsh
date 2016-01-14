#!/bin/tcsh
source ./code/code.main/custom-tcshrc   # shell settings

##
## USAGE: chipseq-heatmaps.tcsh OUTPUT-DIR PARAMETER-SCRIPT MATRICES-BRANCH [OBJECTS]
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

# filter out inputs
if ($include_input == 'false') set objects = `echo $objects | tr ' ' '\n' | grep -vi input`

# obtain info from sample sheet
set matrices = `./code/read-sample-sheet.tcsh $sheet "$objects" $group yes | awk -v pref=$branch '{print pref"/"$1"/matrix.tsv"}'`
./code/read-sample-sheet.tcsh $sheet "$objects" $group yes | awk -v pref=$branch '{print $2":"$1"\t"pref"/"$1"/matrix.tsv"}' >! $outdir/dataset.tsv

# filter reference loci
cat $matrices | cut -f1 | sort -u | grep -viE "$chrom_excluded" | grep -v '^$' >! $outdir/loci.txt

# create heatmaps
scripts-heatclustering.r -v -o $outdir --row-filter=$outdir/loci.txt $heatmap_params $outdir/dataset.tsv 

# cleanup
#rm -f $outdir/*.RData

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


