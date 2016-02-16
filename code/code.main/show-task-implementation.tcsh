#!/bin/tcsh
#$ -S /bin/tcsh
source ./code/code.main/custom-tcshrc     # customized shell environment

##
## USAGE: show-task-implementation.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

grep ^Rscript __*/run | sed 's/.run:.*results\/.op//' | tr ' ' '\t' | tr -d '"'
echo

foreach d (__*)
  echo `cd $d/inpdirs; ls -1 . | grep -v matrix-estimated`
end

foreach d (__*)
  set d1 = `echo $d | cut -d'-' -f2-`
  grep ^Rscript $d/run | sed 's/\$cmd/.\/code\/hicseq-'"$d1.tcsh/" | sed 's/\$op/'"$d1/g" | sed 's/\$results/results/' | sed 's/\$inpdirs/inpdirs\/*/' 
end
echo


