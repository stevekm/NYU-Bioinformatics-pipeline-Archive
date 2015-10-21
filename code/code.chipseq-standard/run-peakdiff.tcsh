#!/bin/tcsh

##
## USAGE: run-peakdiff.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# set path
set path = (./code/code $path)

# setup
set sheet = inputs/sample-sheet.tsv
set comparisons = inputs/comparisons.tsv
if (! -e $comparisons) then
  scripts-send2err "Warning: comparison file \'$comparisons\' not found."
  exit
endif

# run easydiff
scripts-send2err "=== Running peakdiff ============="
scripts-create-path results/
set threads = 1
set jid = ()
set X = `cat $comparisons | cut -f1`
set Y = `cat $comparisons | cut -f2`
set inpdir = matrices/results
set D = `cd $inpdir; find . -name 'matrix.tsv' | sed 's/\/matrix.[^/]*\/matrix.tsv//' | sort -u | grep 'nbins_1/' | sed 's/^\.\///'`
foreach params (params/params.*.tcsh)
  set params_name = `echo $params | sed 's/^params\/params\.//' | sed 's/\.tcsh$//'`
  set m = 1
  while ($m <= $#X)
    set x = $X[$m]
    set y = $Y[$m]
    set comparison = ${x}.${y}
    if ( (`cat $sheet | grep -v '^#' | cut -f4 | tr '|' '\n' | grep -c "^$x"'$'` == 0) || (`cat $sheet | grep -v '^#' | cut -f4 | tr '|' '\n' | grep -c "^$y"'$'` == 0) ) then
      scripts-send2err "Warning: sample $x or $y not found, skipping..."
    else
      scripts-send2err "-- [$params] comparing $x to $y..."
      foreach d ($D)
        set xfiles = `cat $sheet | grep -v '^#' | cut -f2,4 | tr '|' ' ' | key_expand | cols -t 1 0 | grep "^$x	" | cut -f2 | awk -v d=$inpdir/$d '{print d"/matrix."$0"/matrix.tsv"}'`
        set yfiles = `cat $sheet | grep -v '^#' | cut -f2,4 | tr '|' ' ' | key_expand | cols -t 1 0 | grep "^$y	" | cut -f2 | awk -v d=$inpdir/$d '{print d"/matrix."$0"/matrix.tsv"}'`
        set outdir1 = results/easydiff.$params_name/$d
        scripts-create-path $outdir1/logs/
        set outdir2 = $outdir1/easydiff.$comparison
        set jid = ($jid `scripts-qsub-wrapper $threads ./code/peaks-diff $outdir2 $params "$xfiles" "$yfiles"`)
      end
    endif
    @ m ++
  end
end

scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "$jid"
scripts-send2err "Done."





