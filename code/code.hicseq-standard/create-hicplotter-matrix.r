#!/bin/Rscript

# Rscript create-hicplotter-matrix.r OUT-MATRIX INPUT-MATRIX 

# Get the arguments
args <- commandArgs(trailingOnly = TRUE)

out_mat <- args[1]
input_mat <- read.table(sprintf("%s", args[2]), header=TRUE, row.names=1, check.names=FALSE)
# Write to output
write.table(input_mat, file=sprintf("%s", out_mat), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
