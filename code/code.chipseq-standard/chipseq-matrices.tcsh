#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # shell settings (must be included in all scripts)

##
## USAGE: create-matrix.tcsh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH OBJECT(S)
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
set peak_branch = $branch
source ./code/code.main/scripts-read-job-vars $peak_branch "$objects" "genome genome_dir branch"
set aln_branch = $branch
set branch = $peak_branch

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- window = $win"
scripts-send2err "- flank = $flank"
scripts-send2err "- nbins = $nbins"

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine input files
set alignments = `echo $samples | tr ' ' '\n' | awk -v d=$aln_branch '{print d"/"$0"/alignments.bam"}'`
set ref_regions = `ls -1 $branch/*/peaks.bed`

# create ref.bed
scripts-create-path $outdir/
set ref = $outdir/ref.bed
if ($win > 0) then
  set p = `scripts-create-temp $outdir`
  cat $ref_regions | gtools-regions bed >! $p
  cat $genome_dir/genome.bed | gtools-regions win -s $win -d `echo $win/4 | bc` | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' | gtools-overlaps subset -i $p >! $ref
  rm -f $p
else
  cat $ref_regions | gtools-regions bed | scripts-sortbed | gtools-regions link | gtools-regions center | gtools-regions pos -op 1 | gtools-regions shiftp -5p -$flank -3p +$flank | gtools-regions shiftp -5p 0 -3p -1 | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' >! $ref
endif

# check reference labels  
if (`gtools-regions reg $ref | cut -f1 | sort | uniq -d | wc -l` > 0) then
  scripts-send2err "Error: unique labels required in the $ref file."
  exit 1
endif

# create matrix
set threads = $#alignments
set alignments = `echo $alignments | tr ' ' ','`
gtools-threaded matrix -v -i -p $threads --overlap-op hits -nbins $nbins -rpkm -o $outdir/matrix.tsv $alignments $ref

# cleanup
#rm -f $ref

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


