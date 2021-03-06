#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains-armatus.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set inpdir = $branch/$object
set workdir = $outdir/work
mkdir -p $workdir

set est_matrices = `cd $inpdir; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded"`            # TODO: inpdir should not have both tsv and RData, check!
foreach est_mat ($est_matrices)
  scripts-send2err "Processing matrix $est_mat..."
	set chr = `echo $est_mat | cut -d'.' -f2`
	
	# extract matrices from RData file
	if (`echo $est_mat | grep -c '\.RData$'` == 1) then
    Rscript ./code/hic-matrix.r matrices -v -o $workdir/tmp $inpdir/$est_mat
  else
    mkdir -p $workdir/tmp
    cat $inpdir/$est_mat >! $workdir/tmp/matrix.k=001.tsv
  endif
  
  # convert matrices and run armatus            # TODO: parallelize this?
  foreach mat ($workdir/tmp/matrix.*.tsv)
    set pref = `basename $mat .tsv | sed 's/^matrix\.//'`.$chr
    set inpmat = $pref.inpmat.gz
  	sed '1d' $mat | sed 's/:/	/' | sed 's/-/	/' | sed 's/NA/0/g' | gzip >! $workdir/$inpmat
    set p = `pwd`
    cd $workdir
    armatus -i $inpmat $armatus_params -o $pref -m
    rm -f $inpmat
    cd $p
    rm -f $workdir/../Asums.txt
  end
  rm -rf $workdir/tmp
end
	
# Create genome consensus
set kappas = `cd $workdir; ls -1 *.consensus.txt | cut -d'.' -f1 | sort -u` 
foreach k ($kappas)
  cat $workdir/$k.*.consensus.txt | grep -vwE 'chr[MY]' | sed 's/^chr//g' | sed 's/^X/23/g' | sort -k1,1n -k2,2n | sed 's/^23	/X	/g' | sed 's/^/chr/g' >! $workdir/../domains.$k.bed
end

# cleanup
rm -rf $workdir

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


