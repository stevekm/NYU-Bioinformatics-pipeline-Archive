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
#echo "$0 OUTPUT HERE:"
#echo "\n\n"
#
#echo "$argv" " is argv"
#echo "$outdir is outdir"
#echo "$params is params"
#echo "$branch is branch"
#echo "$objects is objects"
#echo "\n\n"

# get the fastq files from the sample sheet
set R1 = `./code/read-sample-sheet.tcsh $sheet $object fastq-r1`
set R2 = `./code/read-sample-sheet.tcsh $sheet $object fastq-r2 | grep -v '^NA$'`
if ("$R1" != "") set R1 = `echo $R1 | tr ',' '\n' | awk -v d=$branch '{print d"/"$0}'`
if ("$R2" != "") set R2 = `echo $R2 | tr ',' '\n' | awk -v d=$branch '{print d"/"$0}'`
set R1 = `echo $R1 | tr ' ' ','`
set R2 = `echo $R2 | tr ' ' ','`

# check file extensions
set ext = `echo $R1 | tr ',' '\n' | sed 's/\.gz$//' | sed 's/.*\.//' | sort -u`
if ($ext != "fastq") then
  scripts-send2err "Error: only fastq files can be used as input files."
  exit 1
else 
  scripts-send2err "Using fastq files..."
endif

echo "$R1 is R1"
echo "$R2 is R2"
echo "\n\n"

# run FastQC on every file found
foreach i ("$R1" "$R2")
  echo "$i" " is the file to be processed"
  
  # set the output dirname for FastQC
  set fastqcdir = `basename "$i" |  sed 's/\.\(.*\)//'`
  echo $fastqcdir " is the fastq dir name"
  
  # full outdir path, need a subdir for each fastq file
  echo "$outdir"/"$fastqcdir" " is the outdir subdir"
  mkdir -p "$outdir"/"$fastqcdir"
  fastqc --threads "$threads" --nogroup --outdir "$outdir"/"$fastqcdir" "$i"
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."

