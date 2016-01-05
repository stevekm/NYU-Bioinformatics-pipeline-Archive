#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compare-boundaries.tcsh OUTDIR PARAMS BRANCH OBJECT1 OBJECT2 
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = $4
set object2 = $5

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$object1 $object2" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Get the object directories
set domaindir1 = $branch/$object1
set domaindir2 = $branch/$object2

# Concatenate blacklisted regions
cat $black_lists | gtools-regions bed >! $outdir/blacklist.bed

# kappas
set kappas = `ls -1 $domaindir1/domains.k=*.bed $domaindir2/domains.k=*.bed | sed 's/.*domains\.//' | sed 's/\.bed$//' | sort | uniq -d`   # TODO: adapt for boundaries.k=*.bed
set pref = domains
echo "SAMPLE1 SAMPLE2 KAPPA N-SAMPLE1 N-SAMPLE2 N-UNION N-COMMON OVERLAP" | tr ' ' '\t' >! $outdir/comparison.tsv
foreach k ($kappas)
  # find common boundaries	
  set input_file1 = $domaindir1/$pref.$k.bed
  set input_file2 = $domaindir2/$pref.$k.bed
  set boundary_file1 = $outdir/boundaries1.$k.bed
  set boundary_file2 = $outdir/boundaries2.$k.bed

  # Create the boundaries for sample1
  ( gtools-regions pos -op 5p $input_file1 ; gtools-regions pos -op 3p $input_file1 ) | gtools-regions shiftp -5p -$flank_dist -3p +$flank_dist | gtools-overlaps subset -inv $outdir/blacklist.bed | scripts-sortbed >! $boundary_file1

  # Create the boundaries for sample2
  ( gtools-regions pos -op 5p $input_file2 ; gtools-regions pos -op 3p $input_file2 ) | gtools-regions shiftp -5p -$flank_dist -3p +$flank_dist | gtools-overlaps subset -inv $outdir/blacklist.bed | scripts-sortbed >! $boundary_file2

  # First, take the union
  cat $boundary_file1 $boundary_file2 | scripts-sortbed | gtools-regions link >! $outdir/union.$k.bed

  # Then, take the intersection
  cat $boundary_file1 | gtools-overlaps intersect $boundary_file2 | scripts-sortbed | gtools-regions link >! $outdir/intersection.$k.bed

  # Save the common boundaries in the outdir
  cat $outdir/union.$k.bed | gtools-overlaps subset $outdir/intersection.$k.bed >! $outdir/common_boundaries.$k.bed

  # Report stats
  set n1 = `cat $boundary_file1 | wc -l`
  set n2 = `cat $boundary_file2 | wc -l`
  set n_union = `cat $outdir/union.$k.bed | wc -l`
  set n_common = `cat $outdir/common_boundaries.$k.bed | wc -l`
  set a = `echo $n_common/$n_union | bc -l`
  echo "$object1 $object2 $k $n1 $n2 $n_union $n_common $a" | tr ' ' '\t' >> $outdir/comparison.tsv

  # Create the Venn diagram
  # Rscript ./code/create-venn-diagram.r $n1 $n2 $n_common "$object1" "$object2" $outdir
end

# cleanup
rm -f $outdir/blacklist.bed

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."





