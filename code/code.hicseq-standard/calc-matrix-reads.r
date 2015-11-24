#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

##
## USAGE: calc-matrix-reads.r MATRIX
##


args <- commandArgs(trailingOnly=T)
if (length(args)==0) {
  cat("USAGE: calc-matrix-reads.r MATRIX-1 ...\n")
  quit(save="no")
}

n = 0
for (f in args) {
  mat = as.matrix(read.table(f,header=T))
  n = n + sum(mat[upper.tri(mat,diag=TRUE)])
}
write(n,stdout())

quit(save='no')


