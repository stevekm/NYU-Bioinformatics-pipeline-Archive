#!/usr/bin/Rscript

# Source the tool
# source("/ifs/home/cl3011/ROTATION_3/Resources/Software/TopDom.R")

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {print("USAGE: Rscript run-topdom.r INPUT-MATRIX WINDOW-SIZE TOMDOM-PATH"); quit(save="no")}

# Give the matrix and the window size as arguments
input_matrix <- args[1]
window_size  <- as.numeric(args[2])
topdompath   <- args[3]

# Source the tool
source(sprintf("%s", topdompath))

kappa <- strsplit(input_matrix, "[.]")[[1]][1]
chrom <- strsplit(input_matrix, "[.]")[[1]][2]
prefix <- paste(kappa, chrom, sep=".")

# Run TopDom with the input matrix
# the window size and the prefix
# as arguments
TopDom(matrix.file=sprintf("%s", input_matrix), window.size=window_size, sprintf("%s", prefix))

