#!/bin/tcsh

##
## USAGE: fetch-src-code.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

unalias cp 

foreach x (core.h core.cpp gzstream.h gzstream.c)
  rsync -au ~/Code/basic/$x .
end

foreach x (gtools-intervals.h gtools-intervals.cpp gtools-regions.cpp gtools-overlaps.cpp gtools-scans.cpp gtools-threaded.cpp gtools-hic.cpp)
  rsync -au ~/Code/genomic-tools/$x .
end

foreach x (key-expand.cpp mergeuniq.cpp matrix.cpp matrix2.cpp cols.cpp rows.cpp table.cpp vectors.cpp)
  rsync -au ~/Code/tools/$x tools-$x
end

