#!/bin/tcsh

##
## USAGE: hicseq-annotate-tables.tcsh OUTPUT-DIR HIC-TABLES GENE-ANNOTATION-REG LOCI-OF-INTEREST-REG
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set hic_tables = ($2)
set annot_reg = $3
set loci_reg = ($4)

# Create directory
if (-e $out) then
  send2err "Warning: directory $out already exists, overwriting..."
else
  mkdir -p $out
endif

# Get program directory
set prog_dir = `readlink -f $0 | sed 's/\/[^/]\+$//'`

# Print parameters
send2err "Program directory = $prog_dir"
send2err "Output directory = $out"
send2err "Hi-C tables = $hic_tables"
send2err "Annotation regions = $annot_reg"
send2err "Loci of interest = $loci_reg"

# Create reg file for all bins found in the tables
cat $hic_tables | grep ^chr | cut -f1,2 | tr '\t' '\n' | sort -u | cols -t 0 0 | sed 's/:/ + /' | sed 's/-/ /' | cols -t 1 0 >! $out/bin.reg

# Merge loci-of-interest regions
send2err "Merging loci-of-interest regions..."
echo -n '' >! $out/loci.reg
foreach l ($loci_reg)
  set lname = `echo $l | sed 's/.*\///' | sed 's/\.[^.]\+$//' | tr '/' '_'`
  cat $l | genomic_regions reg | cut -f2 | sed "s/^/$lname	/" >> $out/loci.reg
end
  
# Annotate bins by loci of interest
cat $out/loci.reg | genomic_overlaps overlap -i -label -t '|' $out/bin.reg | cut -f1 | replace_with_tab '|' | cols -t 1 0 | sort -u | mergeuniq -merge -t , >! $out/bin.loci.tsv

# Annotate bins by gene names
genomic_regions reg $annot_reg | genomic_overlaps overlap -i -label -t '|' $out/bin.reg | cut -f1 | replace_with_tab '|' | cols -t 1 0 | sort -u | mergeuniq -merge -t , >! $out/bin.gene.tsv

# Join annotations into a single table
joint -a1 -o 1.1 1.2 2.2 -e N/A $out/bin.gene.tsv $out/bin.loci.tsv >! $out/bin.annotated.tsv

# Collecting results
send2err "Collecting results..."
echo `cat $hic_tables | head -1` locus1-genes locus1-marks locus2-genes locus2-marks | tr ' ' '\t' >! $out/table.annotated.tsv
cat $hic_tables \
  | grep ^chr \
  | sort -k1,1 | joint -1 1 - $out/bin.annotated.tsv \
  | sort -k2,2 | joint -1 2 - $out/bin.annotated.tsv \
  | sed 's/^\([^\t]*\t\)\([^\t]*\t\)/\2\1/' \
  | sort -u \
>> $out/table.annotated.tsv

# Cleanup
send2err "Done."


