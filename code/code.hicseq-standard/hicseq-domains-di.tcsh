#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains-di.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH SAMPLE
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set sample = $4

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$sample" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set inpdir = $branch/$sample
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
  
  # convert matrices and run di
  foreach mat ($workdir/tmp/matrix.*.tsv)
    echo $mat
    set pref = `basename $mat .tsv | sed 's/^matrix\.//'`.$chr
    set inpmat = $pref.matrix.txt
    ./code/create-di-matrix.tcsh $mat $workdir/$inpmat
  end
  rm -rf $workdir/tmp
end

# Get the different kappas and create dirs
set kappas = `cd $workdir; ls -1 *.matrix.txt | cut -d'.' -f1 | sort -u`
foreach k ($kappas)
  mkdir -p $workdir/$k
  mv $workdir/$k*.matrix.txt $workdir/$k
end 

# Get into each of the kappa directories 
# and calculate the DI
foreach k ($kappas)
  ./code/calculate-di.tcsh $params domains.$k.bed $workdir/$k
end

# cleanup
rm -rf $workdir

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


