#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: create-sample-sheet.tcsh GENOME={hg19,mm10}
##
## FUNCTION: create sample sheet automatically from input files in fastq directory
##

# process command-line inputs
if ($#argv != 1) then
  grep '^##' $0
  exit
endif

set genome = $1

# create sample sheet
set inpdir = fastq
set sheet = sample-sheet.tsv
echo "sample group fastq-r1 fastq-r2 genome enzyme" | tr ' ' '\t' >! $sheet
foreach sample (`cd $inpdir; ls -1d *`)
  scripts-send2err "Importing sample $sample..."
  set group = `echo $sample | cut -d'-' -f-2`
  set fastq1 = `cd $inpdir; ls -1 $sample/*_R1.fastq.gz`
  set fastq2 = `cd $inpdir; ls -1 $sample/*_R2.fastq.gz`
  if (`echo $sample | tr '-' '\n' | grep -ic HindIII` == 1) then
    set enzyme = HindIII
  else if (`echo $sample | tr '-' '\n' | grep -ic NcoI` == 1) then
    set enzyme = NcoI
  else 
    set enzyme = NA
  endif
  echo "$sample\t$group\t`echo $fastq1 | tr ' ' ','`\t`echo $fastq2 | tr ' ' ','`\t$genome\t$enzyme" >> $sheet
end

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

