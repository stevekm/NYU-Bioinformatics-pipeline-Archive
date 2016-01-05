#!/bin/tcsh

##
## USAGE: run-peakdiff.tcsh
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# setup
set sheet = inputs/sample-sheet.tsv
set grp_col = `cat $sheet | head -1 | tr '\t' '\n' | grep -n '^group$' | cut -d':' -f1`
set comparisons = inputs/comparisons.tsv
if (! -e $comparisons) then
  scripts-send2err "Warning: comparison file \'$comparisons\' not found."
  exit 1
endif

# run easydiff
scripts-send2err "=== Running peakdiff ============="
scripts-create-path results/
set threads = 1
set jid = ()
set X = `cat $comparisons | cut -f1`
set Y = `cat $comparisons | cut -f2`
set inpdir = inpdirs/matrices/results
set D = `cd $inpdir; find . -name 'job.sh' | sed 's/\/[^/]\+\/job.sh//' | sort -u | grep 'nbins_1/' | sed 's/^\.\///'`
foreach params (params/params.*.tcsh)
  set params_name = `echo $params | sed 's/^params\/params\.//' | sed 's/\.tcsh$//'`
  set m = 1
  while ($m <= $#X)
    set x = $X[$m]
    set y = $Y[$m]
    set comparison = ${x}.${y}
    if ( (`cat $sheet | sed '1d' | cut -f$grp_col | tr '|' '\n' | grep -c "^$x"'$'` == 0) || (`cat $sheet | sed '1d' | cut -f$grp_col | tr '|' '\n' | grep -c "^$y"'$'` == 0) ) then
      scripts-send2err "Warning: sample $x or $y not found, skipping..."
    else
      scripts-send2err "-- [$params] comparing $x to $y..."
      foreach d ($D)
        set xfiles = `cat $sheet | sed '1d' | cut -f1,$grp_col | tr '|' ' ' | tools-key-expand | tools-cols -t 1 0 | grep "^$x	" | cut -f2 | awk -v d=$inpdir/$d '{print d"/"$0"/matrix.tsv"}'`
        set yfiles = `cat $sheet | sed '1d' | cut -f1,$grp_col | tr '|' ' ' | tools-key-expand | tools-cols -t 1 0 | grep "^$y	" | cut -f2 | awk -v d=$inpdir/$d '{print d"/"$0"/matrix.tsv"}'`
        set outdir1 = results/easydiff.$params_name/$d
        set outdir2 = $outdir1/easydiff.$comparison
        set jid = ($jid `scripts-qsub-wrapper $threads ./code/chipseq-peakdiff.tcsh $outdir2 $params $inpdir/$d "$xfiles" "$yfiles"`)
      end
    endif
    @ m ++
  end
end

scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "$jid"
scripts-send2err "Done."





