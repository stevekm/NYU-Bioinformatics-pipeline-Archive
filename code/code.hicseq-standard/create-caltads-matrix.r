#!/bin/Rscript

# Check for required packages
if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages("reshape2")}

# Load the required library
library(reshape2)

# Get the command line arguments
args <- commandArgs(trailingOnly=TRUE)

out_matrix <- args[1]
input_matrix <- args[2]
chr <- args[3]

# Load the data
m <- read.table(sprintf("%s", args[2]), header=TRUE, row.names=1)
size <- dim(m)[1] - 1 
colnames(m) <- c(0:size)  
rownames(m) <- c(0:size)

m$cols <- colnames(m)

# Create long format
m1 <- melt(m, id.vars="cols")

# Write to output
#out <- paste0(chr,".",chr,"_",chr,".txt")
#out1 <- paste(outdir,out,sep="/")
write.table(m1,file=sprintf("%s",out_matrix),sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)
