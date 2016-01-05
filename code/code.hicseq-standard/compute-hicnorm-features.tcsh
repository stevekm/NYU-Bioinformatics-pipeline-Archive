#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: compute-hicnorm-features
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

set BINSIZES = (5 10 20 40 100 200 1000)   # in kb

foreach enzyme_info (NcoI:CCATGG HindIII:AAGCTT)
  set enzyme = `echo $enzyme_info | cut -d':' -f1`
  set motif = `echo $enzyme_info | cut -d':' -f2`
  echo "Processing $enzyme ($motif)..."

  foreach w ($BINSIZES)
    echo "-- processing w=$w..."
    ( cat genome.bed | grep '^chr[0-9]' | tools-cols -t 0 0 1 2 | sed 's/^chr//' | sort -n | cut -f2- ; cat genome.bed | grep chrX ) | gtools-regions win -s ${w}000 -d ${w}000 | gtools-regions reg | tools-cols -t 1 1 | sed 's/ + /:/' | sed 's/ /-/' >! genome.w=${w}kb,d=${w}kb.reg
    set x = `make_temp_vec 3`
    cat $enzyme.effective.reg | gtools-overlaps coverage -i genome.w=${w}kb,d=${w}kb.reg | sort >! $x[1]
    cat genome.w=${w}kb,d=${w}kb.reg | gtools-overlaps overlap -i -label -t '|' $enzyme.gc.reg | cut -f1 | replace_with_tab '|' | tools-mergeuniq -merge | tools-vectors m -n 6 | sort >! $x[2]
    cat genome.w=${w}kb,d=${w}kb.reg | gtools-overlaps overlap -i -label -t '|' $enzyme.mappability.reg | cut -f1 | replace_with_tab '|' | tools-mergeuniq -merge | tools-vectors m -n 6 | sort >! $x[3]
    joint -a1 -e0 -o 1.1 1.2 2.2 $x[1] $x[2] | joint -a1 -e0 -o 1.1 1.2 1.3 2.2 - $x[3] | sed 's/:/\t/' | sed 's/-/\t/' | sort -k1,1 -k2,2g | sed 's/\t/:/' | sed 's/\t/-/' >! $enzyme.w=${w}kb.features.txt
    rm -f $x
  end
  
end

echo Done.


