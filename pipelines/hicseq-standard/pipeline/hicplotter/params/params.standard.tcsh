#!/bin/tcsh

source ./inputs/params/params.tcsh

# regions to plot
set regions = "chr10:10000000-15000000 chr10:5000000-10000000 chr10:6000000-8000000"
set bedgraphs = "/ifs/home/cl3011/ROTATION_3/Resources/Software/HiCPlotter/data/CUTLL1-DMSO-CTCF-rep1.bedGraph"
set bedgraph_labels = "CUTLL1-DMSO-CTCF"
set highlight = 0
set highlight_bed = "/ifs/home/cl3011/wrappers/hicplotter_test/highlight.bed"
set loop_bed = "/ifs/home/cl3011/wrappers/hicplotter_test/looplist.bed"
set histogram_type = 0     # Provide comma-separated values if you have more than one
set histogram_height = 100 # Provide comma-separated values if you have more than one
set fileheader = 0         # Either 1 or 0 (header / no header)
set insulation_score = 0   # Either 1 or 0 (include insulation index or not)
set hicplotter_path = "/ifs/home/cl3011/wrappers/hicplotter_test"
set hicplotter_matrix_path = "/ifs/home/cl3011/wrappers/hicplotter_test/code2"
set pdftk_path = "/usr/bin"


