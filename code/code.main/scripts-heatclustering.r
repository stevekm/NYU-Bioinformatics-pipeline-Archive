#!/usr/bin/Rscript


# global variables

VERSION = '2.0a'


#
#  FUNC run_gsea
#
run_gsea <- function(outdir,Clist,fdr,release_dir,use_refg) 
{
  if (use_refg==TRUE) { refg = paste(outdir,'/ref.g',sep='') } else { refg = '' }
  n = length(Clist);
  for (i in 1:n) {
    # create gene list
    fout = paste(outdir,'/cluster_',formatC(i,width=3,format='d',flag='0'),'.g',sep='');
    write.table(Clist[[i]],file=fout,row.names=F,col.names=F,quote=F);
    # run GSEA
    cmd = paste('func_gene_set_with_fdr ',fout,' ',fdr,' ',release_dir,' ',refg,sep='');
    write(cmd,stderr())
    system(cmd);
  }
}



#
#  FUNC gsea_plot
#
gsea_plot <- function(f,main='',type='go',ntop=30,maxn=200) 
{
  if (length(readLines(f))<=1) { plot.new(); title(main); return(); }
  m = read.table(f,sep='\t',skip=1,row.names=NULL,comment.char='',quote='',stringsAsFactors=FALSE)
  if (type=='go') { colnames(m) = c('id','n','qval','pval','val','descr')
  } else if (type=='msigdb') { colnames(m) = c('n','qval','pval','val','descr','descr2') }
  m = m[m$n<maxn,,drop=FALSE]
  k = min(nrow(m),ntop)
  if (k==0) { plot.new(); title(main); return(); }
  m = m[1:k,,drop=FALSE]
  val = rep(NA,ntop)
  val[1:k] = rev(-log10(m$pval))
  lab = rep(NA,ntop)
  lab[1:k] = as.character(rev(m$descr))

  par(las=2) # make label text perpendicular to axis
  par(mar=c(5,18,4,2)) # increase y-axis margin.
  barplot(val,main=main,horiz=T,names.arg=lab,cex.names=0.5,col='green4',border=NA)
}



multi_merge <- function(x)
{
  if (length(x)<2) return(x)
  y = x[[1]]
  for (k in 2:length(x)) y = merge(y,x[[k]],by='ID')
  return(y)
}



# ###################################################################################################################################################
#  OPERATION = ANALYZE
# ###################################################################################################################################################
op_heatmap <- function(cmdline_args) 
{
  # usage
  usage = "\
  heatclustering.r [OPTIONS] DATASET-FILE\
\
Function:\
  Create clustered heatmap for genomics data. Use --help for list of options.\
\
Input files:\
  DATASET-FILE                     tab-separated file containing labels (column #1) and matrix filenames (column #2)\
\
Output files:\
  clustering.RData                 R binary file containing data, including intermediate steps\
  clustering.tif                   clustered heatmap image\
  clustering.txt                   gene cluster assignments\
  profiles.pdf                     cluster profiles\
  cluster_*.g                      list of genes in cluster\
  cluster_*.g.out.go.txt           gene ontology enriched terms (ordered by q-value)\
  cluster_*.g.out.go.gene.txt      gene ontology enriched terms per gene\
  cluster_*.g.out.msigdb.txt       MSigDB enriched genesets (ordered by q-value)\
  cluster_*.g.out.msigdb.gene.txt  MSigDB enriched genesets per gene\
"

  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default=\"%default\"]."),
    make_option(c("-P","--palette"), default="", help="File containing palette colors [default=\"%default\"]."),
    make_option(c("--normalize"), default="rowmax", help="Matrix normalization: none, max, mean, rowmax (default), rowrange, or rowmean."),
    make_option(c("--log2"), default="false", help="Log2 transformation: false (default) or true."),
    make_option(c("--row-filter"), default="", help="File containing labels of rows to be used [default=\"%default\"]."),
    make_option(c("--n-best"), default=2000, help="Keep top rows ranked by range [default=%default]."),
    make_option(c("--diff"), default=0.5, help="Filter out rows that have small differences from min to max [default=%default]."),
    make_option(c("--use-short-names"), action="store_true",default=FALSE, help="Use only second part of the name (after the colon)."),
    make_option(c("--nclust"), default=5, help="Number of clusters [default=%default]."),
    make_option(c("--release-dir"), default="", help="ENSEMBL release directory for GSEA sets [default=\"%default\"]."),
    make_option(c("--use-refg"), action="store_true",default=FALSE, help="Use reference set of genes for GSEA analysis [default=FALSE]."),
    make_option(c("--gsea-fdr"), default=0.05, help="False discovery rate cutoff for enrichment analysis [default=%default].")
  );

  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf));
  opt <- arguments$options
  files <- arguments$args
  if (length(files)!=1) { write(usage,stderr()); quit(save='no'); }

  # input arguments
  fdataset = files[1]
  outdir = opt$"output-dir"
  if (opt$verbose) { write("OPTIONS: \n",stderr()); print(opt,file=stderr()); }
  
  # Error checking on input arguments
  if (outdir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(outdir)==FALSE) { dir.create(outdir) } else { write('Warning: output directory already exists, overwriting!',stderr()) }
  if ((opt$palette!='')&(file.exists(opt$palette)==FALSE)) { write('Error: cannot open palette file!',stderr()); quit(save='no'); }
      
  # read sample labels and matrices
  if (opt$verbose) { write('Reading sample labels and matrices...',stderr()) }
  dataset = read.table(fdataset,header=FALSE,sep='\t',check.names=FALSE,row.names=NULL,stringsAsFactors=FALSE)
  colnames(dataset) = c('labels','files')
  if (sum(duplicated(dataset$labels))>0) { write('Error: sample labels must be unique!',stderr()); quit(save='no') }
  Xlabels_long = dataset$labels
  if (opt$"use-short-names"==FALSE) { Xlabels = Xlabels_long } else { Xlabels = as.vector(sapply(Xlabels_long,function(x){strsplit(x,':')[[1]][2]})) }

  n_matrices = length(dataset$files)
  Xraw = lapply(dataset$files,read.table,header=TRUE,sep='\t',check.names=FALSE,row.names=NULL,stringsAsFactors=TRUE)
  n_bins = sapply(Xraw,ncol)-1
  for (k in 1:n_matrices) {
    colnames(Xraw[[k]])[1] = 'ID'
    if (sum(duplicated(Xraw[[k]]$ID))>0) { write(paste('Error: duplicate IDs found in matrix file "',dataset$files[k],'"!',sep=''),stderr()); quit(save='no') }
  }

  # [Yjoined] join matrices
  if (opt$verbose) { write('Joining matrices...',stderr()) }
  Xjoined = multi_merge(Xraw)                       # join all matrices
  IDjoined = as.vector(Xjoined[,1])                 # matrix row IDs (e.g. gene names)
  Yjoined = as.matrix(Xjoined[,2:ncol(Xjoined)])    # joined matrix values-only
  
  # [Ysplit] split joined matrix Y into original groups (will be submatrices of Xraw if row IDs are not identical across all matrices)
  Ysplit = {}
  j = 1; for (k in 1:n_matrices) { jj = j+n_bins[k]-1; Ysplit[[k]] = Yjoined[,j:jj]; j = jj+1 }

  # [Ynorm] normalize joined matrix Y
  if (opt$verbose) { write('Normalizing data...',stderr()); }
  if (opt$normalize=="none") {
    Ynorm = Yjoined
  } else if (opt$normalize=="max") {
    Ynorm = do.call(cbind,lapply(Ysplit,function(b) b/max(b)))
  } else if (opt$normalize=="mean") {
    Ynorm = do.call(cbind,lapply(Ysplit,function(b) b/mean(b)))
  } else if (opt$normalize=="rowmax") {
    Ynorm = do.call(cbind,lapply(Ysplit,function(b) b/apply(b,1,max)))
    Ynorm[is.na(Ynorm)] = 0
  } else if (opt$normalize=="rowrange") {
    Ynorm = do.call(cbind,lapply(Ysplit,function(b) (b-apply(b,1,min))/(apply(b,1,max)-apply(b,1,min))))
    Ynorm[is.na(Ynorm)] = 0
  } else if (opt$normalize=="rowmean") {
    Ynorm = do.call(cbind,lapply(Ysplit,function(b) b/apply(b,1,mean)))
    Ynorm[is.na(Ynorm)] = 0
  }

  # [Ynorm] log2 transform
  if (opt$log2) { Ynorm[Ynorm<=0] = min(Ynorm[Ynorm>0]); Ynorm = log2(Ynorm); Ynorm = Ynorm-min(Ynorm) }

  # [Ymax, Yscore, Ydiff] calculate metrics for each row to be used for filtering
  Ymax = matrix(0,nrow(Ynorm),n_matrices)
  j = 1; for (k in 1:n_matrices) { jj = j+n_bins[k]-1; Ymax[,k] = apply(Ynorm[,j:jj],1,max); j = jj+1 }
  Yscore = apply(Ymax,1,max)-apply(Ymax,1,min)
  Ydiff = 1-apply(Ymax,1,min)/apply(Ymax,1,max)
  
  # [Yfilt] filter by row label, fold-change and range
  Ifiltered = Ydiff>=opt$"diff"
  Ifiltered = Ifiltered&(Yscore>=sort(Yscore,decreasing=TRUE)[min(opt$"n-best",length(Yscore))])
  if (opt$"row-filter"!='') {
    row_labels = read.table(opt$"row-filter",header=FALSE,sep='\t',check.names=FALSE,row.names=NULL,as.is=TRUE,comment.char='')[,1]
    Ifiltered = Ifiltered&(IDjoined %in% row_labels)
  }
  if (opt$verbose) write(paste(sum(Ifiltered),' rows after filtering.',sep=''),stderr())
  Yfilt = Ynorm[Ifiltered,]
  IDfilt = IDjoined[Ifiltered]

  # cluster
  if (opt$nclust>1) {
    # cluster
    if (opt$verbose) { write('Clustering...',stderr()); }
    Ysum = matrix(0,nrow(Yfilt),n_matrices)
    j = 1; for (k in 1:n_matrices) { jj = j+n_bins[k]-1; Ysum[,k] = apply(Yfilt[,j:jj],1,sum); j = jj+1 }
    Cobj = kmeans(Ysum,opt$nclust,iter.max=1000,nstart=10)
    
    # re-order clusters
    if (opt$verbose) { write('Re-ordering clusters...',stderr()); }
#    D2 = as.dist(cor(t(as.matrix(Cobj$centers)))+1)
    D2 = dist(as.matrix(Cobj$centers))
    cluster_order = hclust(D2,method="complete")$order
    Cobj$centers = Cobj$centers[cluster_order,]
    rownames(Cobj$centers) = 1:nrow(Cobj$centers)
    Cobj$cluster = order(cluster_order)[Cobj$cluster]
  } else {
    # sort Yfilt rows by sum of values (TODO: parametrize this?)
    Yfilt_order = order(apply(Yfilt,1,sum),decreasing=TRUE)
    Yfilt = Yfilt[Yfilt_order,]
    IDfilt = IDfilt[Yfilt_order]
    # create a single cluster
    Cobj = {}
    Cobj$cluster = rep(1,nrow(Yfilt))
    Cobj$centers = matrix(apply(Yfilt,2,mean),1)
    rownames(Cobj$centers) = "1"
    colnames(Cobj$centers) = colnames(Yfilt)
    cluster_order = NA
  }

  # create data frame and write clustering results to file
  Cdata = data.frame(ID=IDfilt,cluster=Cobj$cluster)
  write.table(Cdata,file=paste(outdir,'/clustering.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')
  
  # define sample colors
  sample_class = sub(':.*','',Xlabels_long)
  sample_class = factor(sample_class,levels=unique(sample_class))
  n_classes = length(unique(sample_class))
  if (opt$palette=='') {
    samplePalette = rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#00E442", "#0072B2", "#D55E00", "#CC79A7"),10)
  } else {  
    samplePalette = read.table(opt$palette,header=FALSE,sep='\t',check.names=FALSE,row.names=NULL,as.is=TRUE,comment.char='')[,1]
  }
  sampleColors = rep(samplePalette,n_classes/length(samplePalette)+1)
  
  
  # profiles
  if (opt$verbose) { write('Creating profiles...',stderr()); }
  centers = {}        # cluster centers by sample
  Z = t(sapply(1:max(Cobj$cluster),function(c) apply(Yfilt[Cobj$cluster==c,,drop=FALSE],2,mean)))
  j = 1; for (k in 1:n_matrices) { jj = j+n_bins[k]-1; centers[[k]] = Z[,j:jj,drop=FALSE]; colnames(centers[[k]]) = 1:n_bins[k]; j = jj+1 }
  names(centers) = Xlabels
  m = melt(centers)      # data frame
  colnames(m) = c('cluster','bin','val','sample') 
  m$sample = factor(m$sample,levels=names(centers))
  m = merge(m,data.frame(sample=names(centers),class=sample_class))                   # add sample class
  p1 = ggplot(m,aes(x=bin,y=val,group=sample,colour=class))+geom_line()+facet_grid(cluster~.)
  p1 = p1 + 
       scale_colour_manual(values=samplePalette) + 
       theme(axis.text=element_blank(),axis.ticks=element_blank()) + xlab('distance from reference region center') + ylab('normalized read density')
  p2 = ggplot(m,aes(x=bin,y=val,group=sample,colour=class,fill=class))+geom_line()+geom_area()+facet_grid(sample~cluster) 
  p2 = p2 + 
       scale_colour_manual(values=rep('black',n_classes)) + 
       scale_fill_manual(values=sampleColors) + 
       ggtitle('cluster') +
       theme(legend.position='none') + theme(strip.text.y=element_text(size=8,colour='black',angle=0)) +
       theme(axis.text=element_blank(),axis.ticks=element_blank()) +
       xlab('distance from reference region center') + ylab('normalized read density')
  pdf(file=paste(outdir,'/profiles.pdf',sep=''))
  print(p1)
  print(p2)
  dev.off()

  # heatmap
  if (opt$verbose) { write('Creating heatmap...',stderr()); }
  tiff(filename=paste(outdir,'/clustering.tif',sep=''), width=3000, height=3000, compression="lzw", units="px", pointsize=36)  
  layout(matrix(1:(n_matrices+1),1,n_matrices+1),widths=c(0.025,0.975*n_bins/sum(n_bins)))                                 # create layout
  par(mar=c(15,0,0,0),oma=c(2,2,2,2),cex=1)
  a = c(1,sapply(1:opt$nclust,function(x) which(sort(Cobj$cluster,decreasing=TRUE)==x)[1]-1)/(length(Cobj$cluster)-1))     # create cluster color band (rows)
  clustPalette = c('blue','red','green4','cyan','pink','brown','magenta','blue','purple')
  clustColors = rep(clustPalette,length(Cobj$cluster)/length(clustPalette)+1)
  image(matrix(0,1),axes=FALSE,xlim=c(0,1),ylim=c(0,1))
  for (k in 1:opt$nclust) rect(0,a[k+1],1,a[k],col=clustColors[k],border=NA)
  zlim = c(min(Yfilt),max(Yfilt))
  j = 1                                                                                                                    # plot heatmaps, one sample at a time
  for (k in 1:n_matrices) {
    group_palette = colorRampPalette(c('white',sampleColors[sample_class[k]]))(n=100)
    jj = j+n_bins[k]-1
    if (opt$nclust==1) img = t(apply(Yfilt,2,rev)[,j:jj])                  # keep existing order if nclust=1
    else img = t(Yfilt[order(Cobj$cluster,decreasing=TRUE),j:jj])          # otherwise, order by cluster number   # TODO: check this
    image(img,zlim=zlim,xaxt='n',yaxt='n',col=group_palette)
    mtext(text=Xlabels[k],side=1,las=2,cex=1.0)
    j = jj+1
  }
  dev.off();

  # run GSEA analysis
  if (opt$"release-dir"!="") {
    if (opt$verbose) { write('Performing enrichment analysis...',stderr()); }
    Clist = tapply(as.character(Cdata$ID),Cdata$cluster,FUN=list);
    write.table(IDfilt,file=paste(outdir,'/ref.g',sep=''),row.names=F,col.names=F,quote=F,sep='\t');
    run_gsea(outdir,Clist,opt$"gsea-fdr",opt$"release-dir",opt$"use-refg");
    F = dir(outdir,pattern='cluster_[0-9]+.g.out.[^.]+.txt')
    pdf(file=paste(outdir,'/gsea.pdf',sep=''))
    for (f in F) {
      ff = paste(outdir,f,sep='/')
      if (length(grep('go.txt',f))>0) { gsea_plot(ff,main=f,type='go',ntop=30,maxn=200)
      } else if (length(grep('msigdb.txt',f))>0) { gsea_plot(ff,main=f,type='msigdb',ntop=30,maxn=Inf)
      }
    } 
    dev.off()
  }

  # save and exit
  if (opt$verbose) { write('Saving data...',stderr()); }
#  save(Xraw,Xlabels,Ynorm,Yfilt,Ifiltered,IDfilt,Cobj,Cdata,cluster_order,file=paste(outdir,'/clustering.RData',sep=''));
  if (opt$verbose) { write('Done.',stderr()); }
}



  

# ########################
#  MAIN
# ########################

# process command-line arguments
args <- commandArgs(trailingOnly=T)

# install packages
for (p in c("optparse","gplots","ggplot2","reshape")) 
  if (!require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,verbose=FALSE)
  }

# run
op_heatmap(args)

quit(save='no')







