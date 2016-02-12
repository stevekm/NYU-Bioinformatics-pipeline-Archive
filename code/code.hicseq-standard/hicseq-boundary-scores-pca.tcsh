#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-boundary-scores-pca.tcsh OUTDIR PARAM-SCRIPT BOUNDARY-SCORES-BRANCH [OBJECTS]
##

if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Check for number of objects
if ($#objects < 2) then
	scripts-send2err "Error: more than one input objects are required."
	exit 1
endif 

# Generate PCA plots
./code/read-sample-sheet.tcsh $sheet "$objects" $group_var yes | awk '{print $2":"$1}' >! $outdir/labels.tsv
set K = `cd $branch/$objects[1]; ls -1 all_scores.k=*.tsv | cut -d'.' -f2`
set methods = `cat $branch/$objects[1]/all_scores.$K[1].tsv | head -1 | cut -f2-`
set m = 2
foreach method ($methods)
  scripts-send2err "Processing $method..."
  foreach k ($K)
    echo -n "" >! $outdir/data.tsv
    foreach object ($objects)
      cat $branch/$object/all_scores.$k.tsv | cut -f1,$m | sed '1d' | grep -vwE "$chrom_excluded" | sed "s/\t/\t$object\t/" >> $outdir/data.tsv
    end
    cat $outdir/data.tsv | tools-table -c -n 6 | sed 's/ *$//' | tr -s ' ' '\t' >! $outdir/matrix.tsv
    scripts-perform-pca.r -v -o $outdir -L $outdir/labels.tsv --show-text --use-short-names --plain $outdir/matrix.tsv
    cp $outdir/report.qnorm.pdf $outdir/pca.$method.$k.pdf
    rm -f $outdir/data.tsv $outdir/matrix.tsv $outdir/report.*.pdf
  end
  @ m ++
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."




