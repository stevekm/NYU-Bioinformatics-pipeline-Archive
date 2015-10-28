#!/usr/bin/Rscript


# global variables

VERSION = '1.0'



normalize_matrix <- function(D)
{
  D_norm <- normalize.quantiles(D);
  dimnames(D_norm) <- dimnames(D);
  return(D_norm);
}



calc_fdr_cutoff <- function(pos,neg,fdr) 
{
  if (fdr<=0) return(Inf);
  pos <- sort(pos);
  neg <- sort(neg);
  kpos <- 1;
  kneg <- 1;
  while ((kpos<=length(pos))&(kneg<=length(neg))) {
    if ((length(neg)-kneg+1)/(length(pos)-kpos+1)<=fdr) { break; }
    if (pos[kpos]<neg[kneg]) { kpos <- kpos+1; }
    else if (pos[kpos]>neg[kneg]) { kneg <- kneg+1; }
    else { kpos <- kpos+1; kneg <- kneg+1; }
  }
  if (kpos>length(pos)) { y <- 1.01*pos[length(pos)] }
  else { y <- pos[kpos]; }
  return(y);
}



calc_fdr_cutoff_with_bins <- function(values,obs_scores,exp_scores,fc_cutoff,fdr,fdr_bin_size)
{
  # values: values (e.g. RPKMs) as a function of which the fdr will be computed
  # obs_scores: observed scores (e.g. fold-changes across samples) 
  # exp_scores: scores expected by chance (e.g. fold-changes within replicates)
  # fc_cutoff: minimum required fold change
  # fdr: false discovery rate cutoff
  # fdr_bin_size: minimum number of instances per bin where FDR will be computed (as a function of value)

  n = length(values)
  ivalues = order(values)
  value_bin = c()
  t_cutoff = c()
  t_score = rep(0,n)
  m = 2                                    # smoothing parameter, number of micro-bins per bin = 2m+1
  mbin_size = fdr_bin_size/(2*m+1)         # micro-bin size
  mbin_starts = seq(1,n,by=mbin_size)
  k = 1
  for (s in mbin_starts) {
    imbin = ivalues[s:min(n,s+mbin_size-1)]
    bin = max(1,s-m*mbin_size):min(n,s+(m+1)*mbin_size-1)                                   # flank microbin by m*mbin_size
    ibin = ivalues[bin]                                                                     # values inside bin
    t_bin_cutoff = max(calc_fdr_cutoff(obs_scores[ibin],exp_scores[ibin],fdr),fc_cutoff)    # cutoff is determined using all values in bin
    t_score[imbin] = obs_scores[imbin]/t_bin_cutoff                                         # score is only updated in micro-bin
    t_cutoff = c(t_cutoff,t_bin_cutoff,t_bin_cutoff)
	  value_bin = c(value_bin,min(values[imbin]),max(values[imbin]))
    k = k + 1
    if (k%%10==0) write(paste('* ',round(100*k/length(mbin_starts),0),'% complete',sep=''),stderr())
  }
  return(list(value_bin=value_bin,t_cutoff=t_cutoff,t_score=t_score))
}



score <- function(x,y,method)      # method = { 'mean', 'paired', 'pairwise' }
{
  if (method=='paired') {
    s = mean(x/y)
  } else if (method=='mean') {
    s = mean(x)/mean(y)
  } else if (method=='pairwise') {
    if (length(y)==0) {
      z = combn(x,2)
#      s = mean(c(z[1,]/z[2,],z[2,]/z[1,]))
      s = mean(z[1,]/z[2,])
    } else {
      z = expand.grid(x,y)
      s = mean(z[,1]/z[,2])
    }
  } else {
    s = NA
  }
  
  return(s)
}


my_ttest <- function(x,y,alternative,paired)
{
  return(0)   # TODO: fix this
  t.test(x,y,alternative='two.sided',paired=paired)$'p.value'
}


diff_peaks.calc <- function(D,signal_cols,ref_cols,fdr,fdr_bin_size,fold_cutoff,method)
{
  write('Initializing...',stderr())
  diffpeaks = {}
  diffpeaks$signal_cols = signal_cols
  diffpeaks$ref_cols = ref_cols
  diffpeaks$fdr = fdr
  diffpeaks$fdr_bin_size = fdr_bin_size
  diffpeaks$fold_cutoff = fold_cutoff

  write('Enforcing positive lower bound on matrix values...',stderr())
  lbound = min(D[D>0])
  D[D<lbound] = lbound 
  diffpeaks$D = D
  
  write('Computing fold-changes between samples...',stderr())
  diffpeaks$gain = apply(D,1,function(x) score(x=x[signal_cols],y=x[ref_cols],method=method))
  diffpeaks$loss = apply(D,1,function(x) score(x=x[ref_cols],y=x[signal_cols],method=method))

  write('Computing p-values...',stderr())
  diffpeaks$pval = apply(log2(D),1,function(z) my_ttest(z[signal_cols],z[ref_cols],alternative='two.sided',paired=ifelse(method=='paired',TRUE,FALSE)))

  write('Computing fold-changes within replicates...',stderr())
  diffpeaks$loss_bg = apply(D,1,function(x) score(x=x[ref_cols],y=NULL,method='pairwise'))
  diffpeaks$gain_bg = apply(D,1,function(x) score(x=x[signal_cols],y=NULL,method='pairwise'))

  write('Computing FDR on peak gains...',stderr())
  diffpeaks$gain_cutoff = calc_fdr_cutoff_with_bins(apply(D[,ref_cols,drop=FALSE],1,mean),diffpeaks$gain,diffpeaks$gain_bg,fold_cutoff,fdr,fdr_bin_size)
  write('Computing FDR on peak losses...',stderr())
  diffpeaks$loss_cutoff = calc_fdr_cutoff_with_bins(apply(D[,signal_cols,drop=FALSE],1,mean),diffpeaks$loss,diffpeaks$loss_bg,fold_cutoff,fdr,fdr_bin_size)
  
  diffpeaks$gain_significant = (diffpeaks$gain>=fold_cutoff)&(diffpeaks$gain_cutoff$t_score>=1)
  diffpeaks$loss_significant = (diffpeaks$loss>=fold_cutoff)&(diffpeaks$loss_cutoff$t_score>=1)
  write(paste('Gain = ',sum(diffpeaks$gain_significant),sep=''),stderr())
  write(paste('Loss = ',sum(diffpeaks$loss_significant),sep=''),stderr())

  return(diffpeaks)
}



diff_peaks.plot <- function(diffpeaks,scale)
{
  f = function(z) { z }
  if (scale=='log2') f = log2
  x = apply(diffpeaks$D[,diffpeaks$ref_cols,drop=FALSE],1,mean)
  y = apply(diffpeaks$D[,diffpeaks$signal_cols,drop=FALSE],1,mean)
  vlim = f(c(min(c(x,y)),max(c(x,y))))
  x_lab = paste(colnames(diffpeaks$D)[diffpeaks$ref_cols[1]],' mean (',scale,')',sep='')
  y_lab = paste(colnames(diffpeaks$D)[diffpeaks$signal_cols[1]],' mean (',scale,')',sep='')
  fclim = log2(c(min(c(diffpeaks$gain,diffpeaks$loss)),max(c(diffpeaks$gain,diffpeaks$loss))))
  # x=ref y=sig/ref
  smoothScatter(f(x),log2(diffpeaks$gain),ylim=fclim,xlab=x_lab,ylab='fold-change (log2)',main='signal vs reference')
  lines(f(diffpeaks$gain_cutoff$value_bin),log2(diffpeaks$gain_cutoff$t_cutoff),col='red')
  # x=sig y=ref/sig
  smoothScatter(f(y),log2(diffpeaks$loss),ylim=fclim,xlab=y_lab,ylab='fold-change (log2)',main='reference vs signal')
  lines(f(diffpeaks$loss_cutoff$value_bin),log2(diffpeaks$loss_cutoff$t_cutoff),col='green')
  # x=ref y=sig
  smoothScatter(f(x),f(y),xlim=vlim,ylim=vlim,xlab=x_lab,ylab=y_lab,main='differential peaks')
  igain = diffpeaks$gain_significant
  points(f(x[igain]),f(y[igain]),pch=18,col='red')
  iloss = diffpeaks$loss_significant
  points(f(x[iloss]),f(y[iloss]),pch=18,col='green')
  # boxplots
  boxplot(log2(diffpeaks$gain),log2(diffpeaks$loss),log2(diffpeaks$gain_bg),log2(diffpeaks$loss_bg))
}



diff_peaks.store <- function(x,y,w_signal,w_ref,diffpeaks,out_prefix)
{
  # compute mean/min for reference and signal samples
  y_filt = y[w_signal&w_ref,]
  ref_mean = apply(y_filt[,colnames(y)=='reference'],1,mean)
  ref_min = apply(y_filt[,colnames(y)=='reference'],1,min)
  sig_mean = apply(y_filt[,colnames(y)=='signal'],1,mean)
  sig_min = apply(y_filt[,colnames(y)=='signal'],1,min)

  # save score data
  score_file <- paste(out_prefix,'.score',sep='')
  scores <- cbind(rownames(y_filt),round(log2(diffpeaks$gain),3),diffpeaks$pval,ref_mean,sig_mean,ref_min,sig_min,y_filt)
  colnames(scores)[1] <- 'locus'
  colnames(scores)[2] <- 'fold-change(log2)';
  colnames(scores)[3] <- 'p-value';
  write.table(scores,score_file,quote=F,row.names=F,col.names=T,sep='\t')
  
  # save gain/loss reg files
  write.table(scores[diffpeaks$gain_significant,c(1,2),drop=FALSE],paste(out_prefix,'.gain',sep=''),quote=F,row.names=F,col.names=F,sep='\t');
  write.table(scores[diffpeaks$loss_significant,c(1,2),drop=FALSE],paste(out_prefix,'.loss',sep=''),quote=F,row.names=F,col.names=F,sep='\t');
  
  # save outlier data
  outlier_file <- paste(out_prefix,'.outliers',sep='');
  d <- rbind(y[!w_signal,],y[!w_ref,]);
  outliers <- round(d,digits=6);
  outliers <- cbind(rownames(d),outliers);
  colnames(outliers)[1] <- 'locus';
  write.table(outliers,outlier_file,quote=F,row.names=F,col.names=T,sep='\t');
  
  # create RData file
  save(x,y,w_signal,w_ref,diffpeaks,file=paste(out_prefix,'.RData',sep=''))
}





remove_outliers <- function(D,outlier_prob,scale)
{
  w <- 1:nrow(D)
  if (ncol(D)==1) return(w)

  # enforce lower bound
  lbound = min(D[D>0])
  D[D<lbound] = lbound 

  # scale
  f = function(z) { z }
  if (scale=='log2') f = log2

  signal_z <- f(D)
  vlim = c(min(signal_z[signal_z>-Inf]),max(signal_z))
  signal_label <- colnames(D)[1]
  smoothScatter(signal_z,xlab=paste(signal_label,' #1 (',scale,')',sep=''),xlim=vlim,ylim=vlim,ylab=paste(signal_label,' #2 (',scale,')',sep=''),main='replicate reproducibility');

  if (outlier_prob>0) {
    n_sample = 20000;
    i_extrema <- as.vector(c(which(signal_z[,1]==min(signal_z[,1]))[1],which(signal_z[,1]==max(signal_z[,1]))[1]));
    i <- c(i_extrema,sample(nrow(signal_z),n_sample,rep=T))
    fit <- loess(signal_z[i,2] ~ signal_z[i,1],span=0.5,degree=1);
    x <- sort(signal_z[i,1]);
    lines(x,predict(fit,x),col='magenta');
    r <- signal_z[,2]-predict(fit,signal_z[,1]);
    w <- dnorm(r,mean(r,na.rm=T),sd(r,na.rm=T))>outlier_prob;
    points(signal_z[!w,],pch=19,col='brown');
  } 

  return(w);
}



# ##############################################
#   op_easydiff
# ##############################################
op_easydiff <- function(cmdline_args) 
{
  # usage
  usage = "\
  easydiff.r [OPTIONS] INPUT-MATRIX\
\
Function:\
  Identifies differences between two samples. Use --help for list of options.\
\
Input files:\
  INPUT-MATRIX          tab-separated input data file (columns are samples), use --help for more details\
\
Output files:\
  diff.RData            RData file containing all relevant data structures used for the analysis \
  diff.gain             gains \
  diff.loss             losses \
  diff.outliers         outliers \
  diff.pdf              scatter plots \
  diff.score            all scores and data \
"

  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default=\"%default\"]."),
    make_option(c("--nref"), default=0, help="Number of reference samples [default=%default]."),
    make_option(c("--normalize"), default="none", help="Matrix normalization: none or normq [default=%default]."),
    make_option(c("--scale"), default="none", help="Scale to be used for plotting: none or log2 [default=%default]."),
    make_option(c("--method"), default="paired", help="Method for fold change computations: paired, mean, pairwise [default=%default]."),
    make_option(c("--outlier-prob"), default=0.0, help="Outlier probability cutoff [default=%default]."),
    make_option(c("--fdr-cutoff"), default=0.05, help="False discovery rate cutoff [default=%default]."),
    make_option(c("--fc-cutoff"), default=1.5, help="Fold-change cutoff [default=%default]."),
    make_option(c("--val-cutoff"), default=-Inf, help="Value cutoff [default=%default].")
  )

  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf));
  opt <- arguments$options
  files <- arguments$args
  if (length(files)!=1) { write(usage,stderr()); quit(save='no'); }

  # process input parameters
  data_file = files[1]
  nref = opt$'nref'
  out_dir = opt$'output-dir'
  normalization = opt$'normalize'
  scale = opt$'scale'
  method = opt$'method'
  outlier_prob = opt$'outlier-prob'
  fdr = opt$'fdr-cutoff'
  fold_cutoff = opt$'fc-cutoff'
  val_cutoff = opt$'val-cutoff'

  # check parameters
  if (nref<=0) { write('Error: number of reference samples must be greater than zero!',stderr()); quit(save='no') }

  # read data
  x = as.matrix(read.table(data_file,check.names=F,header=T,row.names=1,sep='\t'))
  ref_cols = 1:nref
  signal_cols = (nref+1):ncol(x)
  sample_labels = c('reference','signal')

  # check parameters
  if (length(signal_cols)!=length(ref_cols)) write('Warning: number of reference samples not equal to number of signal samples!',stderr())

  # create output directory
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no') }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Warning: output directory already exists, results will be overwritten!',stderr()) }
  out_prefix = paste(out_dir,'/diff',sep='')
  image_file <- paste(out_prefix,'.pdf',sep='')

  # set column labels
  colnames(x)[ref_cols] = sample_labels[1]
  colnames(x)[signal_cols] = sample_labels[2]

  # normalize
  if (normalization == 'normq') { 
    y = normalize_matrix(x)
  } else { 
    y = x 
  }

  # filter based on original values
  if (opt$verbose) write('Filtering out low values...',stderr())
  i_filtered = (apply(x[,signal_cols],1,mean)<val_cutoff)&(apply(x[,ref_cols],1,mean)<val_cutoff)
  y = y[i_filtered==FALSE,]

  # setup pdf
  pdf(image_file,width=5,height=7); 
  par(mfrow=c(3,2),cex=0.5,mar=c(4,4,4,4))

  # remove outliers
  if (opt$verbose&(outlier_prob>0)) write('Removing outliers...',stderr())
  w_ref <- remove_outliers(y[,ref_cols,drop=FALSE],outlier_prob=outlier_prob,scale=scale)
  w_sig <- remove_outliers(y[,signal_cols,drop=FALSE],outlier_prob=outlier_prob,scale=scale)
  z = y[w_sig&w_ref,]
  
  # determine FDR bin size
  fdr_bin_size = min(floor(nrow(z)/10),10000)
  if (opt$verbose) write(paste("FDR bin size = ",fdr_bin_size,sep=''),stderr())

  # find differential peaks
  diffpeaks = diff_peaks.calc(D=z,signal_cols=signal_cols,ref_cols=ref_cols,fdr=fdr,fdr_bin_size=fdr_bin_size,fold_cutoff=fold_cutoff,method=method)

  # plot peak differences
  if (opt$verbose) write('Plotting differences...',stderr())
  diff_peaks.plot(diffpeaks,scale=scale)
  dev.off()

  # store results
  if (opt$verbose) write('Storing results...',stderr())
  diff_peaks.store(x,y,w_sig,w_ref,diffpeaks,out_prefix)

  if (opt$verbose) write('Done.',stderr())
}







# ##################################################################
#   MAIN PROGRAM
# ##################################################################


# process command-line arguments
args <- commandArgs(trailingOnly=T)

# install packages
for (p in c('optparse','preprocessCore','MASS')) 
  if (!require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,verbose=FALSE)
  }

# run
op_easydiff(args)

quit(save='no')





