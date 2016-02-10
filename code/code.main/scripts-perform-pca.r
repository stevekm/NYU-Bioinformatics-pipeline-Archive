#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

##
## USAGE: scripts-perform-pca.r [OPTIONS] MATRIX
##

plotPCA <- function(mat,show_text,use_short_names,plain)
{
  pca = prcomp(t(mat))
  if (plain==FALSE) {
    plot(pca)
    heatmap.2(pca$x,scale='column',key=FALSE,trace='none',margins=c(1,5),Colv=FALSE,dendrogram='row',labCol=NULL,cexRow=0.5)
  }
  names = colnames(mat)
  fac = factor(sapply(names,function(x){strsplit(x,':')[[1]][1]}))
  short_names = as.vector(sapply(names,function(x){strsplit(x,':')[[1]][2]}))
  if (show_text) { if (use_short_names) { labels = short_names } else { labels = names } } else { labels = NULL }
  colours = rep(c(brewer.pal(7,"Set1"),brewer.pal(7,"Set2"),brewer.pal(7,"Set3")),nlevels(fac))[1:nlevels(fac)]

  if (plain==TRUE) { 
    pc_list = c(PC2~PC1)
    f = function(pc) {
      xyplot(
        pc, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1,
        panel=function(x, y, ...) { panel.xyplot(x, y, ...); ltext(x=x, y=y, labels=labels, pos=1, offset=0.8, cex=0.5) },
        aspect = "fill", col=colours, scales=list(x=list(at=NULL),y=list(at=NULL)), xlab=NULL, ylab=NULL
      )
    }
  } else { 
    pc_list = c(PC2~PC1,PC3~PC2) 
    f = function(pc) {
      xyplot(
        pc, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1,
        panel=function(x, y, ...) { panel.xyplot(x, y, ...); ltext(x=x, y=y, labels=labels, pos=1, offset=0.8, cex=0.5) },
        aspect = "fill", col=colours,
        main = draw.key(key=list(rect=list(col=colours),text=list(levels(fac)),rep=FALSE,cex=0.5))
      )
    }
  }

  # generate plots
  lapply(lapply(pc_list,f),print)

  return(pca)
}




################################
##### MAIN   ###################
################################


cmdline_args <- commandArgs(trailingOnly=T);

# install packages
for (p in c("optparse")) 
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }

# process command-line arguments
option_list <- list(
  make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
  make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
  make_option(c("-L","--sample-labels"), default="", help="A file containing sample labels in its first column [default \"%default\"]."),
  make_option(c("--show-text"), action="store_true",default=FALSE, help="Display sample label text on PCA plot."),
  make_option(c("--use-short-names"), action="store_true",default=FALSE, help="Use only second part of the name (after the colon)."),
  make_option(c("--plain"), action="store_true",default=FALSE, help="Create only one page (PC1 vs PC2).")
);
usage = 'perform_pca.r [OPTIONS] MATRIX';
  
# get command line options (if help option encountered print help and exit)
arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
opt <- arguments$options;
files <- arguments$args;
if (length(files)!=1) { write('Error: this operation requires an input matrix!',stderr()); quit(save='no'); }

# input arguments
out_dir = opt$'output-dir'
sample_labels = opt$'sample-labels'
show_text = opt$'show-text'

# create output directory
if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Warning: output directory already exists, it will be overwritten!',stderr()) }

# load data
if (opt$'verbose'==TRUE) write('Loading input matrix...',stderr())
x = as.matrix(read.table(files[1],header=T,row.names=1,check.names=FALSE))
if (sample_labels!="") { colnames(x) = t(read.table(sample_labels,header=F)[,1]) }

# load libraries
library("RColorBrewer")
library("lattice")
library("gplots")

# PCA on raw input matrix
if (opt$'verbose'==TRUE) write('Performing PCA on input matrix...',stderr())
pdf(paste(out_dir,'/report.raw.pdf',sep=''))
x_pca = plotPCA(x,show_text=show_text,use_short_names=opt$"use-short-names",plain=opt$"plain")
dev.off()

# PCA on mean-normalized matrix
if (opt$'verbose'==TRUE) write('Performing PCA on mean-normalized matrix...',stderr())
z = t(t(x)/apply(x,2,mean))
pdf(paste(out_dir,'/report.mnorm.pdf',sep=''))
z_pca = plotPCA(z,show_text=show_text,use_short_names=opt$"use-short-names",plain=opt$"plain")
dev.off()

# PCA on quantile-normalized matrix
if (opt$'verbose'==TRUE) write('Performing PCA on quantile-normalized matrix...',stderr())
suppressMessages(library("preprocessCore"))
y = normalize.quantiles(x)
rownames(y) = rownames(x)
colnames(y) = colnames(x)
pdf(paste(out_dir,'/report.qnorm.pdf',sep=''))
y_pca = plotPCA(y,show_text=show_text,use_short_names=opt$"use-short-names",plain=opt$"plain")
dev.off()

#if (opt$'verbose'==TRUE) write('Saving results...',stderr())
#save(x_pca,y_pca,file=paste(out_dir,'/data.RData',sep=''))

if (opt$'verbose'==TRUE) write('Done.',stderr())
quit(save='no')




