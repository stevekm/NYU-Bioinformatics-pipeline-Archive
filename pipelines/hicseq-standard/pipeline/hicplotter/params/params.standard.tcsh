#!/bin/tcsh

source ./inputs/params/params.tcsh

# HiCplotter path
set hicplotter_path = ./code/HiCPlotter2.py

# create bedgraphs for boundary scores
set bscores_branch = ../boundary-scores/results/boundary-scores.by_sample.standard/`echo $branch | sed 's/.*results\///'`
set cell_type = `echo $objects[1] | cut -d'-' -f1`
set f = $bscores_branch/$objects[1]/all_scores.k=001.tsv
set methods = (intra-max DI ratio)
set bedgraphs = ()
set bedgraph_labels = ($methods)
foreach m ($methods)
  set k = `head -1 $f | tr '\t' '\n' | grep -n "^$m"'$' | cut -d':' -f1`
  cat $f | sed '1d' | cut -f1,$k | sed 's/:/\t/' | sed 's/-/\t/' >! $outdir/bscores.$m.bedGraph
  set bedgraphs = ($bedgraphs $outdir/bscores.$m.bedGraph)
end

# add CTCF ChIP-seq
if (-e inputs/data.external/$cell_type/CTCF.bedGraph) then
  set bedgraphs = ($bedgraphs inputs/data.external/$cell_type/CTCF.bedGraph)
  set bedgraph_labels = ($bedgraph_labels CTCF)
endif

# regions to plot
set regions = "chr8:125000000-133000000"
set tiles = "params/regions.bed"
set tiles_labels = "regions"
set highlight = 1
set highlight_bed = "params/highlight.bed"
set fileheader = 0         # Either 1 or 0 (header / no header)
set insulation_score = 0   # Either 1 or 0 (include insulation index or not)



