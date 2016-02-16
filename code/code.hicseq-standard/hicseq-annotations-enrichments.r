#!/bin/Rscript

# Load the required packages
if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages("reshape2", repos='http://cran.us.r-project.org')}
if("corrplot" %in% rownames(installed.packages()) == FALSE) {install.packages("corrplot", repos='http://cran.us.r-project.org')}

library(reshape2)
library(corrplot)
library(plyr)

########## Extract the interactants ###########
split   <- function(x)(unlist(strsplit(x, "[,]")))
#extract <- function(x){x1 <- split(x[8]); x2 <- split(x[10]); x3 <- expand.grid(x1,x2); return(x3)}
# Function I need
f <- function(r) { x = c(as.vector(r["locus1-marks"])); y = c(as.vector(r["locus2-marks"])); x1 = unlist(strsplit(x,"[,]")); y1 = unlist(strsplit(y,"[,]")); expand.grid(x1,y1)}
###############################################

args <- commandArgs(trailingOnly=TRUE)

# Get the arguments (annotations table)
if (length(args) != 3) {print("plot_annotations.R OUTDIR ANNOTATION-FILE CUTOFF"); quit(save="no")}

# Get the annotation file
outdir     <- args[1]
annot_file <- args[2]
cutoff     <- as.integer(args[3])

# Read it
annot <- read.table(sprintf("%s", annot_file), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

# Get the unique loci
N = length(unique(unlist(annot[,1:2])))
# Print the number
sprintf("Number of unique loci: %s", N)

# Sort based on distance-normalized values
annot1 <- annot[order(annot[5],decreasing=TRUE),]
# Select only the top selected ones
annot2 <- annot1[1:cutoff,,drop=FALSE]

# Test if any marks are present in the top-scoring interactions
if (sum((annot2[,"locus1-marks"]!="N/A")|(annot2[,"locus2-marks"]!="N/A"))==0) {
  write("Warning: no marks found in top-scoring interactions.",stderr())
  quit(save='no')
}

# Get a list of interactants
annot2_l <- apply(annot2,1,f)
# Convert to dataframe
annot2_df <- do.call(rbind, annot2_l)
# Count the interactions
value1 <- rep(1,dim(annot2_df)[1])
annot2_df <- cbind(annot2_df,value1)

# Create the raw count matrix
annot2_m <- acast(annot2_df, Var1~Var2, value.var='value1')

# Now order both rownames and colnames alphabetically
annot2_m1 <- annot2_m[order(rownames(annot2_m)),]
annot2_m2 <- annot2_m1[,order(colnames(annot2_m1))]

# Counts of the top pairs
counts_f <- paste(outdir,"top_counts.tsv",sep="/")
write.table(annot2_m2,file=sprintf("%s", counts_f),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

# Now calculate the enrichment
x1 = unique(annot1[,c("locus1","locus1-marks")])
x2 = unique(annot1[,c("locus2","locus2-marks")])
colnames(x1) = colnames(x2) = c("locus","marks")
frequencies <- count(unlist(strsplit(unique(rbind(x1,x2))[,2],',')))

# First find the frequency entries that are represented in the top interactions selected
frequencies_f <- subset(frequencies, frequencies$x %in% rownames(annot2_m2))
colnames(frequencies_f) <- c("ChIP","Count")
freq_f <- paste(outdir,"counts.tsv",sep="/")
write.table(frequencies_f,file=sprintf("%s", freq_f),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

# Create the frequencies vector
counts <- as.vector(frequencies_f[,2])

# Now calculate the enrichment
enrich <- t((annot2_m2/cutoff)/(counts/N))/(counts/N)
#enrich <- t(annot2_m2/counts)/counts
# Remove rows and columns with NA
row_idx <- grep("N/A", rownames(enrich))
col_idx <- grep("N/A", colnames(enrich))
enrich1 <- enrich[-c(row_idx),-c(col_idx)]
enrich_f <- paste(outdir,"enrich.tsv",sep="/")
write.table(enrich1,file=sprintf("%s", enrich_f),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

# write.table(enr)
out <- paste(outdir,"enrichment.pdf",sep="/")
pdf(sprintf("%s",out), useDingbats=FALSE)
corrplot(enrich1, is.corr=FALSE, order="alphabet", method="circle",
         tl.pos="lt", tl.col="blue", cl.pos="b", tl.srt=45, tl.cex=0.8)
dev.off() 

