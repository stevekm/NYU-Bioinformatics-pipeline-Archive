#!/bin/tcsh

##
## USAGE #1: chipseq-align.tcsh OUTPUT-DIR PARAMETER-SCRIPT FASTQ-R1 [FASTQ-R2]
##
## USAGE #2: chipseq-align.tcsh OUTPUT-DIR PARAMETER-SCRIPT BAM-FILE
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

if ($#argv < 3) then
  grep '^##' $0
  exit
endif

set out_dir = $1
set params = $2
if (`echo $3 | grep -c '\.bam$'` == 1) then
  set BAM = $3
  set R1 = 
  set R2 =
else
  set BAM = 
  set R1 = $3
  set R2 = $4
endif
    
# set parameters
source $params
if (! $?NSLOTS) then
  set threads = 8
else
  set threads = $NSLOTS
endif

# create path
scripts-create-path $out_dir

if ($R1 != '') then
  # run aligner
  scripts-send2err "Aligning reads..."
  if ($R2 == "") then
    set input = "-U $R1"
  else
    set input = "-1 $R1 -2 $R2"
  endif
  $aligner --threads $threads $align_params $input | samtools view -q $mapq -@ $threads -Sb1 - | samtools sort -m 10G -@ $threads - $out_dir/alignments_sorted
else
  scripts-send2err "Sorting alignments..."
  samtools sort -m 10G -@ $threads $BAM $out_dir/alignments_sorted
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
 METRICS_FILE=$out_dir/picard_metrics.txt \
 INPUT=$out_dir/alignments_sorted.bam \
 OUTPUT=$out_dir/alignments.bam

# index
samtools index $out_dir/alignments.bam

# stats
scripts-send2err "Computing statistics..."
if ($R1 != '') then
  set R1_files = `echo $R1 | tr ',' ' '`
  set n_reads = `cat $R1_files | gunzip | grep ^@ | wc -l`
else
  set n_reads = `samtools view $BAM | wc -l`
endif
set n_aligned = `samtools view $out_dir/alignments_sorted.bam | wc -l`
set n_unique = `samtools view $out_dir/alignments.bam | wc -l`
echo "Total reads\t$n_reads" >! $out_dir/stats.tsv
echo "Aligned reads\t$n_aligned" >> $out_dir/stats.tsv
echo "De-duplicated alignments\t$n_unique" >> $out_dir/stats.tsv

# create bigwig
scripts-send2err "Creating bigwig file..."
set chr_sizes = `scripts-create-temp`
set bedgraph = `scripts-create-temp`
cat $release/genome.bed | awk '{print $1,$2,$3,$1}' | genomic_regions n >! $chr_sizes
set scale = `echo 1000000/$n_unique | bc -l`
genomeCoverageBed -ibam $out_dir/alignments.bam -scale $scale -bg -g $chr_sizes >! $bedgraph
bedGraphToBigWig $bedgraph $chr_sizes $out_dir/track.bw
rm -f $chr_sizes $bedgraph

# done
rm -f $out_dir/alignments_sorted.bam
scripts-send2err "Done."

