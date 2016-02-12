#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-domains-topdom.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH SAMPLE
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
  
  # convert matrices and run caltads
  foreach mat ($workdir/tmp/matrix.*.tsv)
    set pref = `basename $mat .tsv | sed 's/^matrix\.//'`.$chr
    set inpmat = $pref.matrix.txt
    Rscript ./code/create-topdom-matrix.r $mat $workdir/$inpmat 
    set tpath = `readlink -f $topdompath`
    set p = `pwd`
    cd $workdir
    Rscript $p/code/run-topdom.r $inpmat $winsize $tpath
    rm -f $inpmat
    cd $p
  end
  rm -rf $workdir/tmp
end
	
# Create genomewide domains
set kappas = `cd $workdir; ls -1 *.domain | cut -d'.' -f1 | sort -u` 
foreach k ($kappas)
  cat $workdir/$k.*.domain | grep -Ev 'chr	' | grep -E "	domain	" | cut -f1,3,5 | sed 's/^chr//' | sed 's/^X	/23	/g' | sort -k1,1n -k2,2n | sed 's/^23	/X	/' | sed 's/^/chr/' | scripts-sortbed >! $workdir/../domains.$k.bed
end

# cleanup
# rm -rf $workdir

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


