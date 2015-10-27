#!/bin/bash

# this script copies the results of a chip-seq analysis over to a Results dir

# call script like this:
# scripts-copy-chipseq.sh {input/project dir full paths} {results dir full path}
# EXMPLE:
# ~/pipeline-master/code/scripts-copy-chipseq.sh ~/projects/Sussel-Giselle-ChIPseq-2015-10-13/ /ifs/data/sequence/results/external/columbia/2015-10-13/chip-seq/

# you can check the arguments here:
# echo $1
# echo $2

# things to copy:
# bigwigs
# peaks; narrow and broad
# PCA
# peaks table

###
# ls -1 --ignore="job*" $1/chipseq-standard/alignments/results/*/* # use this to ignore the job files


### COPY THE BIG WIGS
# make dir
mkdir -p $2/alignment/bigwigs

# find the bigwig files
BIG_WIGS=$(ls $1/chipseq-standard/alignments/results/*/*.bw)
for i in $BIG_WIGS; do
  # get the sample name from its dir name
  q=$(basename $(dirname $i ) );
  
  # strip out things we don't need in the name
  z=${q/*./};
  
  # copy the file with new name (don't overwrite if already present!)
  cp -n $i $2/alignment/bigwigs/$z.bw
done


### COPY THE PEAKS
# make the dirs
mkdir -p $2/peaks/{narrow,broad}

# find the braod peaks
BROAD_PEAKS=$(ls $1/chipseq-standard/peaks/results/peaks.macs_broad/*/peaks.bed)
for i in $BROAD_PEAKS; do
  #echo $i;
  q=$(basename $(dirname $i ) );
  z=${q/*./};
  # echo $z
  cp -n $i $2/peaks/broad/$z.bed
done

# get the broad peak scores
BROAD_PEAKS_SCORES=$(ls $1/chipseq-standard/peaks/results/peaks.macs_broad/*/peak-scores.bed)
for i in $BROAD_PEAKS_SCORES; do
  #echo $i;
  q=$(basename $(dirname $i ) );
  z=${q/*./};
  cp -n $i $2/peaks/broad/${z}_scores.bed
done


# find the narrow peaks
NARROW_PEAKS=$(ls $1/chipseq-standard/peaks/results/peaks.macs_narrow/*/peaks.bed)
for i in $NARROW_PEAKS; do
  #echo $i;
  q=$(basename $(dirname $i ) );
  z=${q/*./};
  # echo $z
  cp -n $i $2/peaks/narrow/$z.bed
done

# get the narrow peak scores
NARROW_PEAKS_SCORES=$(ls $1/chipseq-standard/peaks/results/peaks.macs_narrow/*/peak-scores.bed)
for i in $NARROW_PEAKS_SCORES; do
  #echo $i;
  q=$(basename $(dirname $i ) );
  z=${q/*./};
  cp -n $i $2/peaks/narrow/${z}_scores.bed
done

rsync -a --exclude "logs" $1chipseq-standard/pca/results/ $2pca

