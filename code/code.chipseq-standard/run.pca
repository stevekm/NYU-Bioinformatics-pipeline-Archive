#!/bin/tcsh
#$ -S /bin/tcsh

##
## USAGE: run.pca
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

# set path
set path = (./code/code $path)

# process sample sheet
set sheet = inputs/sample-sheet.tsv                                                                         # TODO: the code below will not work for multi-group assignments
set groups_noinput = `cat $sheet | grep -v '^#' | cut -f4 | tr '|' '\n' | sort | uniq | grep -vwi input`
if ($#groups_noinput <= 1) then
  set awk_filter = '$0'                        # include inputs only if no more than one group (besides inputs) is available
else
  set awk_filter = 'tolower($2)!~/input/'      # remove all inputs from PCA
  scripts-send2err "Warning: Excluding all Inputs from PCA."
endif

# run PCA
scripts-send2err "=== Running PCA ============="
scripts-create-path results/logs/
set inpdir = matrices/results
set D = `cd $inpdir; find . -name 'matrix.*.tsv' | sed 's/\/matrix.*.tsv$//' | sort -u | grep 'nbins=1/' | sed 's/^\.\///'`
set jid = ()
foreach d ($D)
  set out = results/`echo $d | sed "s/matrices\./pca./"`
  scripts-send2err "- out = $out"
  set mat = $out/matrix.tsv
  set lab = $out/labels.txt
  if (! -e $out) then
    scripts-create-path $out/
    if (! -e $mat) then
      cat $sheet | grep -v '^#' | awk $awk_filter | cut -f2,4 | awk '{print $2":"$1}' >! $lab
      set files = `cat $sheet | grep -v '^#' | awk $awk_filter | cut -f2 | awk -v pref=$inpdir/$d/ '{print pref"matrix."$1".tsv"}'`
      scripts-multi-paste $files | grep -v '^chr[MYX]' >! $mat                                                              # TODO: exclude chrM/X/Y?
    else
      scripts-send2err "Warning: $mat already exists, skipping..."
    endif  
    set jid = ($jid `scripts-qsub-run results/logs/job.pca 1 scripts-perform-pca.r -v -o $out -L $lab --show-text --use-short-names $mat`)
  else 
    scripts-send2err "Warning: $out already exists, skipping..."
  endif
end

scripts-send2err "Waiting for all jobs to finish..."
scripts-qsub-wait "$jid"
scripts-send2err "Done."





