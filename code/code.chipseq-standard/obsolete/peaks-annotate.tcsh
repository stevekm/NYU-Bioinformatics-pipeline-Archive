#!/bin/tcsh

##
## USAGE: peaks-annotate.tcsh PARAMS PEAKS REFERENCE
##

# shell settings (must be included in all scripts)
source ./code/code.main/custom-tcshrc

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

set params = $1
set peaks = ($2)
set ref = $3

set proximal = 1000
set upstream = 100000
set downstream = 100000
set annot = `scripts-create-temp`
set table1 = `scripts-create-temp`
set table2 = `scripts-create-temp`

# peak info
( echo "PEAK-ID\tPEAK-LOCUS\tPEAK-SIZE"; \
  cat $peaks | genomic_regions bed | awk '{ print $4"\t"$1":"$2"-"$3"\t"$3-$2 }' | sort \
) >! $table1

# proximal territories
( cat $ref | genomic_regions reg | sed 's/^/genebody:/'; \
  cat $ref | genomic_regions pos -op 5p | genomic_regions shiftp -5p -$proximal -3p +$proximal | genomic_regions reg | sed 's/^/proximal:/'; \
  cat $ref | genomic_regions pos -op 5p | genomic_regions shiftp -5p -$upstream -3p -$proximal | genomic_regions reg | sed 's/^/upstream:/'; \
  cat $ref | genomic_regions pos -op 5p | genomic_regions shiftp -5p +$proximal -3p +$downstream | genomic_regions reg | sed 's/^/downstream:/' \
) >! $annot

# peaks, genes and offsets
( echo "PEAK-ID\tGENE-ID\tREGION\tDISTANCE"; \
  cat $peaks | genomic_regions reg |  genomic_overlaps offset -i -label -c $annot | tr ':' '\t' | awk '{ print $3"\t"$2"\t"$1"\t"$4 }' | sort \
) >! $table2

# join tables
scripts-join-tables $table1 $table2 

# cleanup
rm -f $annot $table1 $table2




