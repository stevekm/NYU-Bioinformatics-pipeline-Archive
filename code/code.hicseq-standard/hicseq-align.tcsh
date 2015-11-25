#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # customize shell environment

##
## USAGE: hicseq-align.tcsh OUTPUT-DIR PARAM-SCRIPT FASTQ-DIR SAMPLE
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set fastq_dir = $3
set sample = $4

# set parameters
send2err "Setting parameters..."
source $params
if (! $?NSLOTS) then
  set threads = 16
else
  set threads = $NSLOTS
endif

# determine input fastq filenames
set fastq1 = `cat $sheet | awk -v s=$sample '$1==s' | cut -f3 | tr ',' '\n' | awk -v d=$fastq_dir '{print d"/"$0}'`
set fastq2 = `cat $sheet | awk -v s=$sample '$1==s' | cut -f4 | tr ',' '\n' | awk -v d=$fastq_dir '{print d"/"$0}'`

# create path
scripts-create-path $out

# align
if ($aligner == 'gtools') then                    ## Aligner = gtools

  if ( ($#fastq1 > 1) || ($#fastq2 > 1)) then
    scripts-send2err "Error: gtools_hic align does not allow multiple read pair files."
    exit
  endif
  gtools_hic align -v --work-dir $out/tmp -p $threads $align_params --bowtie-index $genome_index $fastq1 $fastq2 | samtools view -T $genome_index.fa -b1 - >! $out/alignments.bam
  rm -rf $out/tmp


else if ($aligner == 'bowtie2') then              ## Aligner = bowtie2

  bowtie2 -p $threads $align_params -x $genome_index --sam-nohead -U `echo $fastq1 | tr ' ' ','` -S $out/alignments.R1.sam
  bowtie2 -p $threads $align_params -x $genome_index --sam-nohead -U `echo $fastq2 | tr ' ' ','` -S $out/alignments.R2.sam
  paste -d'\n' $out/alignments.R1.sam $out/alignments.R2.sam | samtools view -T $genome_index.fa -b1 - >! $out/alignments.bam
  rm -f $out/alignments.R1.sam $out/alignments.R2.sam


else
  scripts-send2err "Error: unknown aligner $aligner."
  exit
endif


# save variables
set >! $out/job.vars.tsv




