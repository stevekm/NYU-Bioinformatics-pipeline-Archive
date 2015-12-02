#!/bin/tcsh

source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains-armatus.tcsh OUTPUT-DIR INPUT-DIR GAMMA
##

# Load required module
module load armatus/2014-05-19

set outdir = $1
set input = $2
set gamma = $3

set workdir = $outdir/work
mkdir -p $workdir

set matrices = `cd $input; ls -1 matrix.*.tsv`

foreach matrix ($matrices)
	# Convert the matrix
	set base = `basename $matrix .tsv`
	sed '1d' $input/$matrix | sed 's/:/	/' | sed 's/-/	/' | gzip >! $workdir/$base.gz
	set chr = `echo $base | sed 's/matrix.//'`
	# Go to the output directory
	set p = `pwd`
	cd $workdir
	armatus -i $base.gz -g $gamma -o $chr -m
	rm -f $base.gz
	cd $p
end
	
# Create genome consensus
cat $workdir/*.consensus.txt | grep -viE 'chrM|chrY' | sed 's/^chr//g' | sed 's/^X/23/g' | sort -k1,1n -k2,2n | sed 's/^23	/X	/g' | sed 's/^/chr/g' >! $workdir/../domains.bed

