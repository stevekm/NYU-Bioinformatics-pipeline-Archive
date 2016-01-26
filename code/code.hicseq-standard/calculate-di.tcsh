#!/bin/tcsh

##
## USAGE: calculate_DI.tcsh PARAMS-FILE OUTPUT-DOMAINS INPUT-SAMPLE
##

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

# Get the sample
set params = $1
set domains = $2
set sample = $3

# Get the parameters
source $params

# Get into the sample directory
set p = `pwd`
cd $sample

# Create symbolic link to IC
ln -s $domaincallpath DI_scripts

#Calculate the directionality index
foreach matrix (`ls -1 *matrix.txt`)
  perl DI_scripts/perl_scripts/DI_from_matrix.pl $matrix $bin_size $window_size $faipath >>! DI_temp.txt
end

# Sort the DI file properly
# in order to create the input for HMM
cat DI_temp.txt | grep -vwE 'M|Y' | sed 's/^X	/23	/g' | sort -k1,1n -k2,2n >! DI.txt
rm -rf DI_temp.txt

# Create symbolic links for what is required for the MATLAB run
ln -s $domaincallpath/HMM_calls_test.m
ln -s $domaincallpath/required_modules

# Run MATLAB with DI.txt as input
echo 'Running MATLAB'
nice matlab < HMM_calls_test.m > dumpfile
echo 'MATLAB part complete'

# Get 7 column format
perl DI_scripts/perl_scripts/file_ends_cleaner.pl HMM_outfile.txt DI.txt | perl DI_scripts/perl_scripts/converter_7col.pl >! hmm_7colfile

# Create hmm_7colfile for each chromosome
set chroms = `cut -f1 hmm_7colfile | uniq`

# Create hmm for each chromosome
foreach chrom ($chroms)
  cat hmm_7colfile | grep "^$chrom	" >! $chrom.hmm
end

# Find the domains for each chrom
foreach file (`ls -1 *.hmm`)
  set chrom = `echo $file | cut -d'.' -f1`
  perl DI_scripts/perl_scripts/hmm_probablity_correcter.pl $file $min $prob $bin_size | perl DI_scripts/perl_scripts/hmm-state_caller.pl $faipath $chrom | perl DI_scripts/perl_scripts/hmm-state_domains.pl >! $chrom.finaldomains
end 

# Combine all domains and in the right order - get the genomewide .bed file
cat *.finaldomains | sed 's/^chr//g' | sort -k1,1n -k2,2n | sed 's/^23	/X	/g' | sed 's/^/chr/g' >! ../../$domains

# Get one step up
cd $p 


