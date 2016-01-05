#!/bin/Rscript

# Check for required packages
# and install
for (package in c("plyr","reshape2","RColorBrewer","corrplot")) {
        if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
}

# Load the required packages
library(plyr)
library(reshape2)
library(corrplot)

args <- commandArgs(trailingOnly = TRUE)


# read the data
data <- data.frame(read.table(sprintf("%s", args[1]), header=FALSE, sep="\t"))
colnames(data) <- c("sample1","sample2","pair","method","lambda","correlation")

# Get the name of directory
dir <- dirname(args[1]) 

# Summarise the data
d1 <- ddply(data, c("sample1", "sample2", "pair", "lambda"), summarise, mean=mean(correlation), min=min(correlation), max=max(correlation), N = length(correlation), sd = sd(correlation), se = sd / sqrt(N))
# Get the summary table
out1 <- paste(dir,"summary.tsv",sep="/")
write.table(d1, file=out1, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
        
# Find unique lambdas
lambdas <- unique(d1$lambda)

out2 <- paste(dir,"correlograms.pdf",sep="/")
pdf(out2)
for (i in 1:length(lambdas)) {
	d2 <- subset(d1, d1$lambda==lambdas[i])
  m <- acast(d2, sample1~sample2, value.var='mean')
  col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
	title <- paste("lambda",lambdas[i],sep="=")
  corrplot(m, method="color", col=col(200), type="full", title=sprintf("%s", title), order="original", addCoef.col = "black", tl.cex=0.75, tl.col="black", tl.srt=45, mar=c(0,0,1,0), diag=FALSE)
}
dev.off()

