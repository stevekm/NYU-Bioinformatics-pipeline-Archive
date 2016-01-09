#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # customize shell environment

##
## USAGE: chipseq-align.tcsh OUTPUT-DIR PARAMETER-SCRIPT FASTQ-BRANCH SAMPLE
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# set parameters
source $params
if (! $?NSLOTS) then
  set threads = 8
else
  set threads = $NSLOTS
endif

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set genome_dir = inputs/genomes/$genome

set BAM = 
set R1 = `./code/read-sample-sheet.tcsh $sheet $object fastq-r1`
set R2 = `./code/read-sample-sheet.tcsh $sheet $object fastq-r2 | grep -v '^NA$'`
if ("$R1" != "") set R1 = `echo $R1 | tr ',' '\n' | awk -v d=$branch '{print d"/"$0}'`
if ("$R2" != "") set R2 = `echo $R2 | tr ',' '\n' | awk -v d=$branch '{print d"/"$0}'`
set R1 = `echo $R1 | tr ' ' ','`
set R2 = `echo $R2 | tr ' ' ','`

# check file extensions
set ext = `echo $R1 | tr ',' '\n' | sed 's/\.gz$//' | sed 's/.*\.//' | sort -u`
if ($#ext > 1) then
  scripts-send2err "Error: multiple input file types found."
  exit 1
else if (($ext != "bam") && ($ext != "fastq")) then
  scripts-send2err "Error: only fastq and bam files can be used as input files."
  exit 1
else if ($ext == "bam") then
  scripts-send2err "Using bam files..."
  set BAM = $R1
  set R1 = 
  set R2 =
else 
  scripts-send2err "Using fastq files..."
endif

# align
if ("$R1" != '') then
  # run aligner
  scripts-send2err "Aligning reads..."
  if ("$R2" == "") then
    set input = "-U $R1"
  else
    set input = "-1 $R1 -2 $R2"
  endif
  $aligner --threads $threads $align_params $input | samtools view -q $mapq -@ $threads -Sb1 - | samtools sort -m 10G -@ $threads - $outdir/alignments_sorted
else
  scripts-send2err "Sorting alignments..."
  samtools sort -m 10G -@ $threads $BAM $outdir/alignments_sorted
endif

# remove duplicate alignments
scripts-send2err "Removing duplicates..."
set picard_root = /local/apps/picard-tools/1.88
java -Xms8G -Xmx16G -jar $picard_root/MarkDuplicates.jar \
 VERBOSITY=WARNING \
 QUIET=true \
 VALIDATION_STRINGENCY=LENIENT \
 MAX_RECORDS_IN_RAM=2500000 \
 REMOVE_DUPLICATES=true \
 ASSUME_SORTED=false \
 CREATE_INDEX=false \
 METRICS_FILE=$outdir/picard_metrics.txt \
 INPUT=$outdir/alignments_sorted.bam \
 OUTPUT=$outdir/alignments.bam

# index
samtools index $outdir/alignments.bam

# stats
scripts-send2err "Computing statistics..."
if ("$R1" != '') then
  set R1_files = `echo $R1 | tr ',' ' '`
  set n_reads = `cat $R1_files | gunzip | grep ^@ | wc -l`
else
  set n_reads = `samtools view $BAM | wc -l`
endif
set n_aligned = `samtools view $outdir/alignments_sorted.bam | wc -l`
set n_unique = `samtools view $outdir/alignments.bam | wc -l`
echo "Total reads\t$n_reads" >! $outdir/stats.tsv
echo "Aligned reads\t$n_aligned" >> $outdir/stats.tsv
echo "De-duplicated alignments\t$n_unique" >> $outdir/stats.tsv

# create bigwig
scripts-send2err "Creating bigwig file..."
set chr_sizes = `scripts-create-temp`
set bedgraph = `scripts-create-temp`
cat $genome_dir/genome.bed | awk '{print $1,$2,$3,$1}' | gtools-regions n >! $chr_sizes
set scale = `echo 1000000/$n_unique | bc -l`
genomeCoverageBed -ibam $outdir/alignments.bam -scale $scale -bg -g $chr_sizes >! $bedgraph
bedGraphToBigWig $bedgraph $chr_sizes $outdir/track.bw
rm -f $chr_sizes $bedgraph

# cleanup
rm -f $outdir/alignments_sorted.bam

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."

