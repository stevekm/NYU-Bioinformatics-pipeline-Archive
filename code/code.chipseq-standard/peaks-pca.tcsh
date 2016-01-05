#!/bin/tcsh

##
## USAGE: peaks-pca OUTPUT-DIRECTORY PEAK-FILES
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

# process command-line inputs
if ($#argv != 2) then
  grep '^##' $0
  exit
endif

set out = $1
set peaks = ($2)

scripts-create-path $out/

# create common reference regions
cat $peaks | scripts-sortbed | gtools-regions link | gtools-regions reg | cut -f2 | tools-rows -number -pref PEAK_ >! $out/ref.reg

# generate matrix
gtools-threaded matrix -v -i -p 1 -nbins 1 --overlap-op value `echo $peaks | tr ' ' ','` $out/ref.reg >! $out/matrix.tsv

# run PCA
head -1 $out/matrix.tsv | tr '\t' '\n' | scripts-skipn 1 | cut -d'/' -f1 | cut -d'.' -f2 | tr '-' ':' | sed 's/:/-/' | sed 's/:/-/' >! $out/matrix.labels.txt
scripts-perform-pca.r -v -o $out -L $out/matrix.labels.txt --show-text $out/matrix.tsv


