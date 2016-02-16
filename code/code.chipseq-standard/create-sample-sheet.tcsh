#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: create-sample-sheet.tcsh GENOME={hg19,mm10} [FRAGMENTATION-LENGTH=auto]
##
## FUNCTION: create sample sheet automatically from input files in fastq directory
##

# process command-line inputs
if ($#argv == 0) then
  grep '^##' $0
  exit
endif

set genome = $1
set fraglen = $2

# create sample sheet
set inpdir = fastq
set sheet = sample-sheet.tsv
echo "sample control group fastq-r1 fastq-r2 genome fragmentation-size" | tr ' ' '\t' >! $sheet

# determine sample names
set samples = `cd $inpdir; ls -1d */*.fastq.gz */*.bam | sed 's/.[^/]\+$//' | sort -u` 

# import sample information
foreach sample ($samples)
  scripts-send2err "Importing sample $sample..."
  set control = NA                                    # this should be filled in manually
  set group = `echo $sample | cut -d'-' -f-3`
  set fastq = `cd $inpdir; ls -1 $sample/*.fastq.gz`
  if ($#fastq == 0) then
    set fastq1 = `cd $inpdir; ls -1 $sample/*.bam`             # use bam files if no fastq.gz files are found
    set fastq2 = NA
  else
    set fastq1 = `echo $fastq | tr ' ' '\n' | grep -E '_R1.fastq.gz|_R1_[0-9]\+.fastq.gz'`
    set fastq2 = `echo $fastq | tr ' ' '\n' | grep -E '_R2.fastq.gz|_R2_[0-9]\+.fastq.gz'`       # TODO: check if fastq1 matches fastq2
    if ($#fastq1 == 0) then
      set fastq1 = ($fastq)       # if no R1/R2 patterns are found, assume it is all R1
      set fastq2 = NA
    endif
  endif
  if ("$fastq2" == "") set fastq2 = NA
  if (($fraglen != "") || ($fraglen != "auto")) then
    set frag = $fraglen
  else
    if (`echo $sample | tr '-' '\n' | grep -iEc '^H2AZ|^H[234]K[0-9]'` == 1) then
      set frag = 150   # histone
    else 
      set frag = 400   # TF
    endif
  endif
  echo "$sample\t$control\t$group\t`echo $fastq1 | tr ' ' ','`\t`echo $fastq2 | tr ' ' ','`\t$genome\t$frag" >> $sheet
end

echo
echo "Your sample sheet has been created! Here is how it looks:"
echo
cat $sheet
echo

echo "Diagnostics: "
echo "Field #1: "
cat $sheet | scripts-skipn 1 | cut -f1 | cut -d'-' -f1 | sort | uniq -c
echo "Field #2: "
cat $sheet | scripts-skipn 1 | cut -f1 | cut -d'-' -f2 | sort | uniq -c
echo "Field #3: "
cat $sheet | scripts-skipn 1 | cut -f1 | cut -d'-' -f3 | sort | uniq -c

