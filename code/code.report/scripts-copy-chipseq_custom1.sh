#!/bin/bash
shopt -s extglob  #Enables extglob
# this script copies the results of a chip-seq analysis over to a Results dir

# call script like this:
# scripts-copy-chipseq.sh {input/project dir full paths} {results dir full path}
# EXMPLE:
# ~/pipeline-master/code/scripts-copy-chipseq_custom1.sh ~/projects/ColumbiaU-Jinsook-2015-10-22/ /ifs/data/sequence/results/external/columbia/2015-10-22/chip-seq/

# you can check the arguments here:
# echo $1
# echo $2

# things to copy:
# bigwigs
# peaks; narrow and broad
# PCA
# peaks table

### COPY THE BIG WIGS
# make dir
mkdir -p $2alignment/bigwigs

# find the bigwig files
BIG_WIGS=$(ls ${1}*/alignments/results/*/*.bw)
for i in $BIG_WIGS; do
  # get the sample name from its dir name
  q=$(basename $(dirname $i ) );
  
  # strip out things we don't need in the name
  z=${q/*./};
  
  # copy the file with new name (don't overwrite if already present!)
  cp -an $i $2alignment/bigwigs/$z.bw
done




### COPY THE PEAKS
# make the dirs
mkdir -p $2peaks/{narrow,broad}

# all the files from the peak results that we will send to the client
ALL_PEAKS_FILES=$(ls ${1}*/peaks/results/peaks.macs_*/*/!(*@(.err|.id|.out|.sh)))
for i in $ALL_PEAKS_FILES; do
  # set broad/narrow based on the peak dir
  if [[ $i == *macs_narrow* ]]; then
    PEAKS_nar_brd="narrow"
  else
  
    if [[ $i == *macs_broad* ]]; then
      PEAKS_nar_brd="broad"
    else
      :
    fi
    
  fi
  
  # get the base dir name
  q=$(basename $(dirname $i ) );
  
  # fix the dir name
  z=${q/*./};
  
  # use the name as part of the filename
  beautiful_file_name=${z}_$(basename $i)

  # copy the file to the proper dir with the new name
  cp -an $i $2/peaks/$PEAKS_nar_brd/$beautiful_file_name

done



# copy the PCA
rsync -a --exclude "logs" --exclude 'job*' --exclude '*win*' ${1}*/pca/results/ $2pca

# copy the Peaks table


if [ -f ${1}*/peaktable ]; then
  rsync -a --exclude "job*" --exclude "ref*" ${1}*/peaktable/results/peaktable.standard/ $2peaktable
else
  :
fi
# allow group users to read and write
chmod g+rw -R $2*
chgrp -Rf results $2* # doesn't work yet..? Fail silently


shopt -u extglob # turn off the extglob