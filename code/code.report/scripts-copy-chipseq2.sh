#!/bin/bash
# USAGE: scripts-copy-chipseq.sh /path/to/project/dir /path/to/output/dir
# this script copies the results of a chip-seq analysis over to a Results dir for a client. 

# call script like this:
# 
# EXMPLE:
# ~/pipeline-master/code/scripts-copy-chipseq2.sh ~/projects/ColumbiaU-Jinsook-2015-10-22 /ifs/data/sequence/results/external/columbia/2015-10-22/chip-seq

shopt -s extglob  #Enables extglob

PROJ_DIR=${1}
OUT_DIR=${2}

# make sure the OUT_DIR exists
mkdir -p $OUT_DIR

# for legacy pipeline support check to see if dir chipseq_analysis exists
# # if exists, add to PROJ_DIR path
if [[ -d $PROJ_DIR/chipseq-standard ]]; then
  PROJ_DIR=${PROJ_DIR}/chipseq-standard
fi

#
##
### COPY THE BIG WIGS
# echo "ls ${1}/*/peaks/results"
# find ${1}/*/alignments/results -type f -name "*.bw"
# ./alignments/results/Nkx2_2-2/track.bw
# check if the $PROJ_DIR has peaks results dir
if [[ -d $PROJ_DIR/alignments/results ]]; then
  # make outdir dir
  mkdir -p $OUT_DIR/alignment/bigwigs
  
  # find the bigwig files
  BIG_WIGS=$(find $PROJ_DIR/alignments/results -type f -name "*.bw")
  for i in $BIG_WIGS; do
    # get the sample name from its dir name
    q=$(basename $(dirname $i ) );

    # strip out things we don't need in the name
    z=${q/*./};

    # copy the file with new name (don't overwrite if already present!)
    cp -an $i $OUT_DIR/alignment/bigwigs/$z.bw
  done
# else
  # echo "${1}/*/peaks/results"
  # echo "NOT FOUND"
  # ls ${1}/*/peaks/results
fi
###
##
#



#
##
### COPY THE PEAKS

if [[ -d $PROJ_DIR/peaks/results ]]; then
  # make the outdir
  mkdir -p $OUT_DIR/peaks
  # find the files we want
  FILES=$(find $PROJ_DIR/peaks/results -type f \( ! -name "*.err" ! -name "*.id" ! -name "*.out" ! -name "*.sh" ! -name "*.branch" \) )
  
  for i in $FILES; do
    # use the filepath as the dirname;
    TMP_NAME=$(echo ${i#*results/} | tr / _ )
    echo $TMP_NAME
    cp -an $i $OUT_DIR/peaks/$TMP_NAME
  done
fi
#


#
##
### copy the PCA
if [ -d $PROJ_DIR/pca/results ]; then
  # make the outdir
  mkdir -p $OUT_DIR/pca
  
  FILES=$(find $PROJ_DIR/pca/results -type f \( ! -name "*.win" ! -name "job*" ! -name "obj.*" ! -name "*.sh" ! -name "*.branch" \) )
  for i in $FILES; do
    # use the filepath as the dirname;
    TMP_NAME=$(echo ${i#*pca.standard/} | tr / _ )
    echo $TMP_NAME
    cp -an $i $OUT_DIR/pca/$TMP_NAME
  done
  
  # this preserves the dir structure..
  # rsync -a --exclude "logs" --exclude 'job*' --exclude '*win*' $PROJ_DIR/pca/results/ $OUT_DIR/pca
  
  
fi
# copy the Peaks table
if [ -d $PROJ_DIR/peaktable/results ]; then
  # make the outdir
  mkdir -p $OUT_DIR/peaktable
  
  FILES=$(find $PROJ_DIR/peaktable/results -type f \( ! -name "*.win" ! -name "job*" ! -name "obj.*" ! -name "*.sh" ! -name "*.branch" \) )
  for i in $FILES; do
    # use the filepath as the dirname;
    TMP_NAME=$(echo ${i#*peaktable.standard/} | tr / _ )
    echo $TMP_NAME
    cp -an $i $OUT_DIR/peaktable/$TMP_NAME
  done

  # rsync -a --exclude "job*" --exclude "ref*" $PROJ_DIR/peaktable/results/peaktable.standard/ $OUT_DIR/peaktable
fi
###
##
#

### PERMISSIONS
# allow group users to read and write
chmod g+rw -R $OUT_DIR*
chgrp -Rf results $OUT_DIR* # doesn't work yet..? Fail silently


shopt -u extglob # turn off the extglob