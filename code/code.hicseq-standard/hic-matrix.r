#$ -S /usr/bin/Rscript


# global variables

VERSION = '2.6'


# row correlations (1D)
cor1D <- function(x,y,method) { sapply(seq.int(dim(x)[1]),function(i) cor(as.vector(x[i,]),as.vector(y[i,]),method=method)) }
 
# row correlations (2D)
cor2D <- function(x,y,method) { sapply(seq.int(dim(x)[1]),function(i) cor(as.vector(x[i,,]),as.vector(y[i,,]),method=method)) }


# print_options
print_options <- function(opt)
{
  write("OPTIONS: \n",stderr())
  sink(stderr())
  print(opt)
  sink(stdout())
}

 
# contact matrix stats
matrix_stats <- function(x)
{
  p = 0.99
  q0 = quantile(x,log2(seq(from=2^p,to=2,by=0.001)))
  q1 = quantile(x[upper.tri(x)],log2(seq(from=2^p,to=2,by=0.001)))
  plot(q0,type='l',col='red',lwd=2); 
  lines(q1,col='blue',lwd=2)

}


# compute betas using genlasso
getb_genlasso <- function(y,n)
{
  out = fusedlasso2d(y,verbose=T,maxsteps=maxsteps);
  b = out$beta[,seq(1,ncol(out$beta),length.out=n)];
  return(b);
}

# sample genlasso solution path
sample_fusedlasso2d <- function(y,maxsteps,nsample)
{
  out = fusedlasso2d(y,verbose=T,maxsteps=maxsteps);
  I = seq(1,length(out$lambda),length.out=min(length(out$lambda),nsample));
  Z = list();
  Z[["lambda"]] = out$lambda[I];
  Z[["df"]] = out$df[I];
  Z[["beta"]] = out$beta[,I];
  return(Z);
}


# mark domains on matrix, given boundaries b and estimated contact matrix e
mark_domains <- function(b,e,ntop) 
{
  n = length(b);
  y = matrix(0,n,n);
  z = rollapply(which(b==1),2,by=1,c);   # domains start/stop pairs
  score = rep(0,nrow(z));
  for (i in 1:nrow(z)) { I = z[i,1]:z[i,2]; score[i] = sum(e[I,I]); }
  for (i in order(score,decreasing=T)[1:ntop]) { I = z[i,1]:z[i,2]; y[I,I] = 1; }
  return(y)
}


# count domains given boundaries b
count_domains <- function(b,min_size) 
{
  z = rollapply(which(b==1),2,by=1,c);
  n_domains = 0;
  for (i in 1:nrow(z)) if (z[i,2]-z[i,1]+1>=min_size) n_domains = n_domains + 1;
  return(n_domains)
}



##########################################################################################
#
# FUNCTION: iterative_fused_flsa
# 
# Computes solution path B(gamma,lambda) of column/row-1D fused lasso using flsa.
#
##########################################################################################
iterative_fused_flsa <- function(B,gammas,lambdas,n_iterations,splitCheckSize,verbose)
{
  N = sqrt(length(B[[1]][,1]));
  for (k in 1:length(gammas)) {
    for (j in 1:length(lambdas)) {
      if (verbose) write(paste("Fused1Dalt for gamma=",gammas[k],", lambda=",lambdas[j],"...",sep=''),stderr());
      b = B[[k]][,j];
      for (t in 1:n_iterations) {
        if (verbose) write(paste("- Iteration ",t,"...",sep=''),stderr());
        bb = as.vector(flsa(as.vector(t(matrix(b,N))),lambda1=gammas[k],lambda2=lambdas[j],thr=1e-5,splitCheckSize=splitCheckSize,verbose=verbose));
        e = sum(abs(b-bb));
        write(paste("** error = ",e,sep=''),stderr());
        b = bb;
        if (e<1.0) break;
      }
      B[[k]][,j] = b;
    }
  }
  return(B)
}



#
# CalcPenaltyLossRatio
#
CalcPenaltyLossRatio <- function(y,y_est,lambda,gamma)
{
  if ((lambda==0)&(gamma==0)) return(NA); 
  loss = sum((y-y_est)^2);
  sparsity = sum(abs(y_est));
  fusion_row = sum(abs(y_est[-1,]-y_est[-nrow(y_est),]));
  fusion_col = sum(abs(y_est[,-1]-y_est[,-ncol(y_est)]));
  return((gamma*sparsity + lambda*(fusion_row+fusion_col))/loss);
}



#
# CalcDF
#
CalcDF <- function(X,distance,tolerance)
{
  return(.C("calc_df",X,nrow(X),ncol(X),as.integer(distance),as.double(tolerance),df=integer(1))$df);
}




#
# is_full_matrix
#
is_full_matrix <- function(mat)
{
  return((nrow(mat)==ncol(mat))&&(prod(rownames(mat)==colnames(mat))==1))
}



#
# MatrixImpute
#
MatrixImpute <- function(mat,ignored_rows,ignored_cols) 
{
  if (min(mat,na.rm=TRUE)<0) { write("Error: imputation only works for non-negative matrices!", stderr()); quit(save='no') }
  na_value = -1
  mat[is.na(mat)] = na_value;
  outmat = matrix(.C("impute",as.double(mat),nrow(mat),na_value,y=double(length(mat)))$y,nrow(mat))
  outmat = (outmat+t(outmat))/2
  outmat[ignored_rows,] = 0          # ignored bins are set to zero!
  outmat[,ignored_cols] = 0
  rownames(outmat) = rownames(mat)
  colnames(outmat) = colnames(mat)
  return(outmat)
}



#
# extrema
#
extrema <- function(x,tolerance) 
{
  n = length(x);
  x[x<mean(x)] = min(x);
  #d = x[-1] - x[-n];
  d = (x[-1]-x[-n])/max(abs(x[-1]),abs(x[-n]));
  d[is.nan(d)] = 0;
  d[abs(d)<tolerance] = 0;
  e = c(0,sign(d));
  ee = rep(0,n);
  ee[1] = 1;
  for (k in 2:n) { if ((e[k-1]==1)&(e[k]==-1)) { ee[k-1] = 1 } }
  ee[n] = 1;
  return(ee);
}





#
# MatrixConnectivityDiff (similar to Ren's method?)
#
MatrixConnectivityDiff <- function(x,d)
{
  n = nrow(x);
  c = rep(0,n);
  for (i in (d+1):(n-d)) {
    Iright = seq(i+1,by=1,length.out=d);
    Ileft = seq(i-1,by=-1,length.out=d);
    c[i] = abs(sum(x[i,Iright])-sum(x[i,Ileft]));
  }
  return(c);
}



#
# MatrixBoundaryScores
#
MatrixBoundaryScores <- function(x,distance,d2,skip,ignore)
{
  # initialize
  x[x<0] = 0                           # negative values are set to zero
  x[abs(row(x)-col(x))>distance] = 0   # values beyond distance cutoff are set to zero
  if (d2<=0) d2 = distance
  n = nrow(x)
  labels = c('intra-left','intra-right','intra-max','intra-min','inter','diff','DI','ratio','diffratio','product-max','product-min')
  scores = matrix(0,n,length(labels))
  colnames(scores) = labels
  rownames(scores) = rownames(x)
  
  # compute scores
  for (i in 1:n) {
    # left-of-boundary scores
    Ileft = (i-distance+1):i
    if (min(Ileft)<1) {
      intra_left = intra_left1 = 0
    } else {
      xleft = x[Ileft,Ileft]                                      # left-of-boundary square
      subset = (row(xleft)-col(xleft)>=skip)&(col(xleft)<=d2)     # submatrix
      intra_left = mean(xleft[subset])                            # mean of values in left triangle
      intra_left1 = mean(xleft[1,(skip+1):distance])              # use only 1D, not full triangle
    }
    # right-of-boundary scores
    Iright = i:(i+distance-1)
    if (max(Iright)>n) {
      intra_right = intra_right1 = 0
    } else {
      xright = x[Iright,Iright]                                   # right-of-boundary square
      subset = (row(xright)-col(xright)>=skip)&(col(xright)<=d2)  # submatrix
      intra_right = mean(xright[subset])                          # mean of values in right triangle
      intra_right1 = mean(xright[1,(skip+1):distance])            # use only 1-dimensional vector, not full triangle 
    }
    # scores at boundary
    Ibound = (i-distance+1):(i-1)
    if ((min(Ibound)<1)||(max(Ibound+distance)>n)) {
      inter = 0
    } else {
      inter = mean(x[Ibound,Ibound+distance])*2                   # square at boundary (half of it is zero as it is beyond distance cutoff)
    }
    # left/right difference
    diff = intra_right - intra_left
    diff1 = intra_right1 - intra_left1
    intra_max = max(intra_left,intra_right)                       # max of left/right values
    intra_min = min(intra_left,intra_right)                       # min of left/right values
    # scores
    scores[i,'intra-left'] = intra_left
    scores[i,'intra-right'] = intra_right
    scores[i,'intra-max'] = intra_max                             # max of left/right values
    scores[i,'intra-min'] = intra_min                             # min of left/right values
    scores[i,'inter'] = inter                                     # mean of values at boundary square
    scores[i,'diff'] = diff
    scores[i,'DI'] = ifelse(diff1==0,0,diff1*abs(diff1)/(intra_left1+intra_right1))     # directionality index (Ren 2012)
  }

  # compute regularized ratios
  pseudo = mean(scores[,'inter'])
  scores[,'ratio'] = scores[,'intra-max']/(scores[,'inter']+pseudo)
  scores[,'diffratio'] = scores[,'diff']/(scores[,'inter']+pseudo)

  # product-max
  scores[,'product-max'] = scores[,'intra-max']/max(scores[,'intra-max'])*(1-scores[,'inter']/max(scores[,'inter']))
  scores[,'product-max'] = scores[,'product-max']/max(scores[,'product-max'])      # normalize to [0,1]

  # product-min
  scores[,'product-min'] = scores[,'intra-min']/max(scores[,'intra-min'])*(1-scores[,'inter']/max(scores[,'inter']))
  scores[,'product-min'] = scores[,'product-min']/max(scores[,'product-min'])      # normalize to [0,1]

  return(scores)
}



#
# MatrixDiag
#
MatrixDiag <- function(x,dist)
{
  xx = cbind(x,matrix(0,nrow(x),dist))
  y = t(sapply(1:nrow(x),function(r) xx[r,r:(r+dist)]))
  return(y)
}


#
# MatrixInvDiag
#
MatrixInvDiag <- function(x)
{
  n = nrow(x)
  dist = ncol(x)-1
  y = matrix(0,n,n)                                                # initialize matrix
  for (r in 1:n) { u=min(n,r+dist); y[r,r:u] = x[r,1:(u-r+1)] }    # convert back to original
  y = y+t(y)                                                       # make matrix symmetric
  diag(y) = diag(y)/2
  return(y)
}



#
# MatrixRotate45
#
MatrixRotate45 <- function(x,dist)
{
  if (nrow(x)!=ncol(x)) { write('Error: MatrixRotate45 can only be used with square matrices!',stderr()); quit(save='no'); }
  N = nrow(x)*(dist+1)
  y = .C("rotate45",as.double(x),nrow(x),ncol(x),as.integer(dist),y=double(N))$y
  y = matrix(y,nrow=nrow(x))
  rownames(y) = rownames(x)
  return(y)
}




#
# MatrixInverseRotate45
#
MatrixInverseRotate45 <- function(x)
{
  n = nrow(x)
  dist = ncol(x)
  y = .C("inverse_rotate45",as.double(x),as.integer(n),as.integer(dist),y=double(n*n))$y
  y = matrix(y,n,n)     # convert to matrix
  y = y+t(y)            # make matrix symmetric
  diag(y) = diag(y)/2
  rownames(y) = colnames(y) = rownames(x)
  return(y)
}




#
# calc_counts_per_distance
#
calc_counts_per_distance <- function(mat) 
{ 
  if (is_full_matrix(mat)) {
    n = nrow(mat)
    y = .C("calc_counts_per_distance",as.double(mat),as.integer(n),y=double(n))$y
  } else {
    i = 1; j = 3; while (j<=ncol(mat)) { mat[i,j:ncol(mat)] = NA; i = i + 1; j = j + 2 }
    i = nrow(mat); j = 2; while (j<=ncol(mat)) { mat[i,j:ncol(mat)] = NA; i = i - 1; j = j + 2 }
    y = apply(mat,2,mean,na.rm=TRUE)
  }
  
  return(y)
}




# #####################
#  compare_matrices   #
# #####################

compare_matrices <- function(matrices1,matrices2,distances,method,verbose)
{
  # check solution matrix compatibility
  if (prod(dim(matrices2)==dim(matrices2))!=1) { write("Error: solution matrices have incompatible dimensions!", stderr()); quit(save='no') }
  
  # log2-preprocessing
  if (method=='pearsonlog2') {
    matrices1[matrices1<=0] = min(matrices1[matrices1>0])
    matrices1 = log2(matrices1)
    matrices2[matrices2<=0] = min(matrices2[matrices2>0])
    matrices2 = log2(matrices2)
    method = 'pearson'
  }
  
  # initialize
  n_matrices = dim(matrices1)[1]
  x = matrices1[1,,]
  n_rows = nrow(x)
  n_cols = ncol(x)
  full_matrix = is_full_matrix(x)
  C = matrix(0,nrow=length(distances),ncol=n_matrices)
  rownames(C) = distances

  # compare  
  if (full_matrix==TRUE) {
    distances = unique(sort(c(0,distances)))
    I = col(matrices1[1,,])-row(matrices1[1,,])
    for (i in 2:length(distances)) {
      if (verbose==TRUE) { write(paste('Distance = ',distances[i],':',sep=''),stderr()); }
      J = (I>=distances[i-1])&(I<=distances[i])
      for (k in 1:n_matrices) {
        if (verbose==TRUE) { write(paste('Comparing matrices #',k,'...',sep=''),stderr()); }
        C[i-1,k] = cor(as.vector(matrices1[k,,][J]),as.vector(matrices2[k,,][J]),method=method)
      }
    }
  } else {
    distances[distances>n_cols] = n_cols
    distances = unique(sort(c(0,distances)))
    for (i in 2:length(distances)) {
      if (verbose==TRUE) { write(paste('Distance = ',distances[i],':',sep=''),stderr()); }
      J = distances[i-1]:distances[i]     # select columns according to distance
      for (k in 1:n_matrices) {
        if (verbose==TRUE) { write(paste('Comparing matrices #',k,'...',sep=''),stderr()); }
        C[i-1,k] = cor(as.vector(matrices1[k,,J]),as.vector(matrices2[k,,J]),method=method)
      }
    }
  }
    
  return(C)
}



#
# plot_matrix
#
plot_matrix <- function(mat,main,xlab,ylab,q,qlab,legend_pos,ylim=c(0,1))
{
  n = nrow(mat)
  colors = c('blue','green4','red','magenta','yellow','cyan','brown','pink')
  colors = rep(colors,ceiling(n/length(colors)))[1:n]
  plot(mat[1,],type='l',col=colors[1],ylim=ylim,main=main,xlab=xlab,ylab=ylab,xaxt='n');
  axis(1,at=q,qlab);
  k = 2
  while (k <= n) {
    lines(mat[k,],col=colors[(k-1)%%length(colors)+1])
    k = k + 1
  }
  legend(legend_pos,rownames(mat),lty=rep(1,n),col=colors[1:n],text.col="black",bg="gray90")
}




#
# PreprocessMatrix
#
PreprocessMatrix <- function(x,preprocess,pseudo=1,cutoff=0)       # TODO: add option for dist-norm: custom stats vectors
{
  if (preprocess=="none") {
    y = x
  } else if (preprocess=="max") {
    y = x/max(x)
  } else if (preprocess=="mean") {
    y = x/mean(x)
  } else if (preprocess=="log2") {
    y = log2(x+pseudo)
  } else if (preprocess=="log2mean") {
    y = log2(x+pseudo)
    y = y/mean(y)
  } else if (preprocess=="rank") {
    n = nrow(x)
    y = matrix(rank(x),n)/n^2
  } else if (preprocess=="dist") {
    # TODO: allow custom stats file as input!
    full_matrix = is_full_matrix(x)
    if (full_matrix) {
      counts = calc_counts_per_distance(x)
      D = abs(row(x)-col(x))+1
      expected = matrix(counts[D],nrow(D))
      y = x/expected
    } else {
      expected = calc_counts_per_distance(x)
      y = t(apply(x,1,'/',expected))
    }
    y[is.na(y)] = 0
    y[x<cutoff] = 0
  } else if (preprocess=="distlog2") {
    y = PreprocessMatrix(x,"dist",pseudo=pseudo,cutoff=cutoff)
    y[y<1] = 1
    y = log2(y)
  } else if (preprocess=="zscore") {
    # TODO: allow custom stats file as input!
    full_matrix = is_full_matrix(x)
    if (full_matrix) {
      write('Error: not implemented for full matrices yet!',stderr())
      quit(save='no')
    } else {
      mat = log2(x+pseudo)
      i = 1; j = 3; while (j<=ncol(x)) { mat[i,j:ncol(mat)] = NA; i = i + 1; j = j + 2 }
      i = nrow(mat); j = 2; while (j<=ncol(mat)) { mat[i,j:ncol(mat)] = NA; i = i - 1; j = j + 2 }
      mat_mean = apply(mat,2,mean,na.rm=TRUE)
      mat_sd = apply(mat,2,sd,na.rm=TRUE)
      y = t(apply(t(mat)-mat_mean,2,'/',mat_sd))
    }
    y[is.na(y)] = 0
    y[x<cutoff] = 0
  } else {
    write(paste('Error: unknown preprocessing method "',preprocess,'"!',sep=''),stderr());
    quit(save='no');
  }
  return(y)
}



# ########################
#  OPERATION = ESTIMATE
# ########################

op_estimate <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("--flsa-verbose"), action="store_true",default=FALSE, help="Print FLSA messages."),
    make_option(c("--row-labels"), action="store_true",default=FALSE, help="Input matrix has row labels"),
    make_option(c("-o","--output-file"), default="", help="Output RData file (required) [default=\"%default\"]."),
    make_option(c("-C","--chrom-feature-file"), default="", help="Chromosome feature file (not required) [default=\"%default\"]."),
    make_option(c("--impute"), action="store_true",default=FALSE, help="Impute matrix (only for symmetric matrices; default=replace NAs with zeros)."),
    make_option(c("--ignored-loci"), default="", help="Ignored row and column names found in this list (not required) [default=\"%default\"]."),
    make_option(c("--pseudo"), default=1, help="Pseudocount value to be added to input matrix elements [default=%default]."),
    make_option(c("--min-score"), default=1, help="Minimum HiC score in input matrix [default=%default]."),
    make_option(c("--skip-distance"), default=0, help="Elements near diagonal up to this distance will be set to zero [default=%default]."),
    make_option(c("--preprocess"), default="none", help="Matrix preprocessing: none (default), max, mean, log2, log2mean, rank, dist, distlog2"),
    make_option(c("--min-lambda"), default=0.0, help="Minimum value for parameter lambda [default=%default]."),
    make_option(c("--max-lambda"), default=Inf, help="Maximum value for parameter lambda [default=%default]."),
    make_option(c("--n-lambda"), default=1, help="Number of lambdas [default=%default]."),
    make_option(c("--log2-lambda"), action="store_true",default=FALSE, help="Use log2 scale for lambda range."),
    make_option(c("--gamma"), default=0.0, help="Value for sparsity parameter gamma [default=%default]."),
    make_option(c("--algorithm"), default="fused2D_flsa", help="Algorithm to be used: fused2D_flsa (default) or fused1D_flsa"),
    make_option(c("--zone-size"), default=0, help="Maximum distance from diagonal (number of bins) [default=%default]."),
    make_option(c("--symm"), action="store_true",default=FALSE, help="Make output matrices symmetric (only for fused1D_flsa)."),
    make_option(c("--threshold"), default=1e-05, help="The error threshold used in the algorithm [default=%default]."),
    make_option(c("--split-check-size"), default=1e+09, help="FLSA parameter specifying from which size on, groups of variables are not being checked for breaking up; can be used to reduce computation time; may lead to inaccurate results [default=%default].")
  );
  usage = 'hic-matrix.r estimate [OPTIONS] MATRIX-FILE';

  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf));
  opt <- arguments$options
  if (opt$verbose) print_options(opt) 
  files <- arguments$args;
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }
  
  # input arguments
  fmatrix = files[1]
  fchrom = opt$"chrom-feature-file"
  fignored = opt$"ignored-loci"
  row_labels = opt$"row-labels"
  fout = opt$"output-file"
  pseudo = opt$pseudo
  min_hicscore = opt$'min-score'
  skip_dist = opt$'skip-distance'
  preprocess = opt$preprocess
  min_lambda = opt$"min-lambda"
  max_lambda = opt$"max-lambda"
  n_lambda = opt$"n-lambda"
  log2_lambda = opt$"log2-lambda"
  gammas = opt$gamma
  err_threshold = opt$"threshold"
  splitCheckSize = opt$'split-check-size'
  opt$algorithm = tolower(opt$algorithm)
  algorithm = opt$algorithm
  zone_size = opt$"zone-size"
  rotate45 = zone_size>0
  
  # Error checking on input arguments
  if ((algorithm!='fused2d_flsa')&(algorithm!='fused1d_flsa')) { write(paste('Error: unknown estimation algorithm "',algorithm,'"!',sep=''),stderr()); quit(save='no') }
  if (fout=="") { write('Error: please specify output file!',stderr()); quit(save='no'); }
  
  # read matrix
  if (opt$verbose) { write('Loading matrix...',stderr()); }
  if (row_labels==FALSE) { x = as.matrix(read.table(fmatrix,check.names=F)) } else { x = as.matrix(read.table(fmatrix,row.names=1,check.names=F)) }
  
  # read ignored loci information
  ignored_rows = integer(0)
  ignored_cols = integer(0)
  if (fignored!="") {
    ignored_rows = sort(which(!is.na(match(rownames(x),as.vector(t(read.table(fignored,check.names=F)))))))
    ignored_cols = sort(which(!is.na(match(colnames(x),as.vector(t(read.table(fignored,check.names=F)))))))
  }

  # impute symmetric matrix or replace NAs with zeros
  if ((is_full_matrix(x)==TRUE)&(opt$"impute"==TRUE)) {
    if (opt$verbose) write('Imputing matrix...',stderr())
    x = MatrixImpute(x,ignored_rows,ignored_cols)
    if (opt$verbose) write(paste('Imputed matrix symmetricity = ',sum(x==t(x))/length(x),sep=''),stderr())
  } else {
    x[is.na(x)] = 0
  }
  
  # restrict input matrices to specified distance from diagonal (only applicable to distance-restricted matrices; full matrices are not modified)
  if ((zone_size<=0)||(zone_size>nrow(x))) zone_size = nrow(x)
  if (is_full_matrix(x)==FALSE) x = x[,1:min(ncol(x),zone_size+1)]
 
  # read features
  c = NULL
  if (fchrom!="") { 
    c = read.table(fchrom);
    if (nrow(c)!=nrow(x)) { write('Error: number of rows in contact matrix and chromosome feature matrix must be the same!,stderr()'); quit(save='no'); }
    nreads = sum(x)/1000000;  # in million reads
    wlen = sqrt(c[,4]/1000);  # in kb
    a = 100;  # multiplier
    x = t(x/wlen)/wlen/nreads*a;
    x[is.na(x)] = 0;
  }

  # set ignored rows to zero
  if (opt$verbose) { write(paste('Number of rows set to zero in input matrix (ignored loci) = ',length(ignored_rows),sep=''),stderr()); }
  if (opt$verbose) { write(paste('Number of columns set to zero in input matrix (ignored loci) = ',length(ignored_cols),sep=''),stderr()); }
  x[ignored_rows,] = 0
  x[,ignored_cols] = 0

  # set elements near diagonal to zero
  if (skip_dist>0) {
    D = abs(row(x)-col(x))
    x[D<=skip_dist] = 0
  }
    
  # preprocess input matrix
  y = PreprocessMatrix(x,preprocess=preprocess,pseudo=pseudo,cutoff=min_hicscore)
  
  # create lambda values (zero always included if log2 scale)
  if (max_lambda==Inf) { lambdas = NULL; 
  } else if (log2_lambda==FALSE) { lambdas = unique(seq(min_lambda,max_lambda,length.out=n_lambda));
  } else { lambdas = unique(c(0,2^seq(log2(max(min_lambda,0.01)),log2(max_lambda),length.out=n_lambda-1))); }
  n_lambda = length(lambdas)
  n_matrices = n_lambda
  
  # ##########################################################
  # estimation algorithm = fused/FLSA
  # ##########################################################
  if ((algorithm=='fused2d_flsa')|(algorithm=='fused1d_flsa')) {

    # rotate, if --zone-size>0
    if (rotate45==TRUE) {
      if (opt$verbose) { write('Rotating input matrix...',stderr()); }
      if (opt$verbose) {
        coverage = sum(x[abs(row(x)-col(x))<=zone_size])/sum(x)
        write(paste('coverage = ',coverage,sep=''),file=stderr())
      }
      x = MatrixRotate45(x,zone_size)       # replace input matrix with rotated version
      y = MatrixRotate45(y,zone_size)       # replace preprocessed input matrix with rotated version
      ignored_cols = integer(0)
    }
    
    # if 1D, then convert to vector
    if (algorithm=='fused2d_flsa') { z = y } else if (algorithm=='fused1d_flsa') { z = as.vector(y) }
    
    # estimate
    if (opt$verbose) { write('Estimating contact matrix (rotate=TRUE)...',stderr()); }
    solObj = flsa(z,connListObj=NULL,lambda1=gammas,lambda2=lambdas,verbose=opt$'flsa-verbose',thr=err_threshold,splitCheckSize=splitCheckSize)

    # convert 1D back to 2D if max-lambda!=Inf
    if ((algorithm=='fused1d_flsa')&&(is.null(lambdas)==FALSE)) {
      solObj1 = array(0,dim=c(n_matrices,dim(y)))
      for (i in 1:n_matrices) solObj1[i,,] = matrix(solObj[i,],nrow(y))
      solObj = solObj1
    }
    
  } else {
    write(paste('Error: unknown algorithm "',algorithm,'"!',sep=''),stderr())
    quit(save='no')
  }

  # Post-processing: negative values and ignored rows/columns are set to zero  (TODO: keeping negative values should be an option)
  if (is.null(lambdas)==FALSE) {
    if (opt$verbose) { write('Post-processing estimated matrices...',stderr()); }
    for (i in 1:n_matrices) { 
      solObj[i,ignored_rows,] = 0
      solObj[i,,ignored_cols] = 0
      solObj[i,,][solObj[i,,]<0] = 0
    }
  }

  # save and exit
  if (opt$verbose) { write('Saving data to RData file...',stderr()); }
  save(x,y,opt,ignored_rows,ignored_cols,lambdas,gammas,solObj,file=fout)
  if (opt$verbose) { write('Done.',stderr()) }
  quit(save='no')
}




# ########################
#  GetSolution
# ########################

GetSolution = function(solObj, n_rows, invrotate, lambda, gamma) 
{
  s = flsaGetSolution(solObj,lambda1=gamma,lambda2=lambda)
  if (length(dim(s))==2) { s = s[1,] } else if (length(dim(s))==3) { s = s[1,1,] }
  mat = matrix(s,n_rows)
  if (invrotate==TRUE) return(MatrixInverseRotate45(mat))
  return(mat)
}



# ########################
#  OPERATION = EXTRACT
# ########################

op_extract = function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-file"), default="", help="Output RData file (required) [default=\"%default\"]."),
    make_option(c("--min-lambda"), default=0.0, help="Minimum value for parameter lambda [default=%default]."),
    make_option(c("--max-lambda"), default=Inf, help="Maximum value for parameter lambda [default=%default]."),
    make_option(c("--n-lambda"), default=1, help="Number of lambdas [default=%default]."),
    make_option(c("--log2-lambda"), action="store_true",default=FALSE, help="Use log2 scale for lambda range."),
    make_option(c("--gamma"), default=0.0, help="Value for sparsity parameter gamma [default=%default].")
  )
  usage = 'hic-matrix.r extract [OPTIONS] ESTIMATED-RDATA-FILE';

  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf))
  files <- arguments$args
  new_opt = arguments$options
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }
  
  # read data
  f = files[1]
  if (new_opt$verbose) { write('Loading matrix...',stderr()); }
  load(f)
  if (is.null(opt$lambdas)==FALSE) { write('Error: input RData file should be estimated using max-lambda=Inf.',stderr()); quit(save='no') }
  algorithm = tolower(opt$algorithm)

  # input arguments
  fout = opt$"output-file" = new_opt$"output-file"
  min_lambda = opt$"min-lambda" = new_opt$"min-lambda"
  max_lambda = opt$"max-lambda" = new_opt$"max-lambda"
  n_lambda = opt$"n-lambda" = new_opt$"n-lambda"
  log2_lambda = opt$"log2-lambda" = new_opt$"log2-lambda"
  gammas = opt$gamma = new_opt$gamma
  if (opt$verbose) print_options(opt) 
  
  # Error checking on input arguments
  if (fout=="") { write('Error: please specify output file!',stderr()); quit(save='no'); }
  
  # create lambda values (zero always included if log2 scale)
  if (max_lambda==Inf) { lambdas = NULL; 
  } else if (log2_lambda==FALSE) { lambdas = unique(seq(min_lambda,max_lambda,length.out=n_lambda));
  } else { lambdas = unique(c(0,2^seq(log2(max(min_lambda,0.01)),log2(max_lambda),length.out=n_lambda-1))); }
  n_matrices = n_lambda = length(lambdas)
  n = nrow(y)
  m = ncol(y)

  # create new estimated matrices for specific lambdas/gammas
  if (opt$verbose) { write('Computing solutions for specific lambdas/gammas...',stderr()); }
  solObj0 = solObj
  solObj = array(0,dim=c(n_matrices,n,m))
  for (i in 1:n_matrices) solObj[i,,] = GetSolution(solObj0,n_rows=n,invrotate=FALSE,lambda=lambdas[i],gamma=gammas[1])

  # Post-processing: negative values and ignored rows are set to zero  (TODO: keeping negative values should be an option)
  if (is.null(lambdas)==FALSE) {
    if (opt$verbose) { write('Post-processing estimated matrices...',stderr()); }
    for (i in 1:n_matrices) { solObj[i,ignored_rows,] = 0; solObj[i,,ignored_cols] = 0; solObj[i,,][solObj[i,,]<0] = 0 }
  }

  # save and exit
  if (opt$verbose) { write('Saving data to RData file...',stderr()); }
  save(x,y,z,opt,ignored_rows,ignored_cols,lambdas,gammas,solObj,file=fout)
  if (opt$verbose) { write('Done.',stderr()) }
  quit(save='no')
}

















# ########################
#  OPERATION = PREPROCESS
# ########################

op_preprocess <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("--pseudo"), default=1, help="Pseudocount value to be added to input matrix elements [default=%default]."),
    make_option(c("--row-labels"), action="store_true",default=FALSE, help="Input matrix has row labels"),
    make_option(c("--min-score"), default=1, help="Minimum HiC score in input matrix [default=%default]."),
    make_option(c("--preprocess"), default="none", help="Matrix preprocessing: none (default), max, mean, log2, log2mean, rank, dist, distlog2.")
  )
  usage = 'hic-matrix.r preprocess [OPTIONS] MATRIX';
  
  arguments <- parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf));
  opt <- arguments$options
  if (opt$verbose) print_options(opt) 
  files <- arguments$args;
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }
  
  # input arguments
  fmatrix = files[1]
  row_labels = opt$"row-labels"
  pseudo = opt$pseudo
  min_hicscore = opt$'min-score'
  preprocess = opt$preprocess
  
  # read files
  if (opt$verbose) { write('Loading matrix...',stderr()); }
  if (row_labels==FALSE) { x = as.matrix(read.table(fmatrix,check.names=F)) } else { x = as.matrix(read.table(fmatrix,row.names=1,check.names=F)) }

  # preprocess input matrix
  y = PreprocessMatrix(x,preprocess=preprocess,pseudo=pseudo,cutoff=min_hicscore)

  # print new matrix
  if (opt$verbose) { write("Printing new matrix...",stderr()); }
  write.table(format(y,scientific=TRUE,digits=4),quote=F,sep='\t')

  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no')
}


  



# ##########################
#  OPERATION = NORMALIZE
# ##########################

op_normalize <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-file"), default="", help="Output matrix tsv file (required) [default \"%default\"]."),
    make_option(c("--n-reads"), default=0, help="Total number of (intrachromosomal) reads per sample (required) [default \"%default\"]."),
    make_option(c("--features"), default="", help="Feature file [default \"%default\"]."),
    make_option(c("--ignored-loci"), default="", help="Ignored row and column names found in this list (not required) [default=\"%default\"]."),
    make_option(c("--min-efflen"), default=100, help="Minimum effective length [default \"%default\"]."),
    make_option(c("--min-mappability"), default=0.2, help="Minimum mappability [default \"%default\"]."),
    make_option(c("--scale"), action="store_true",default=FALSE, help="Scale matrix by total number of reads and effective length."),
    make_option(c("--impute"), action="store_true",default=FALSE, help="Impute matrix."),
    make_option(c("--dist"), default=0, help="Distance in number of bins [default \"%default\"].")
  )
  usage = 'hic-matrix.r normalize [OPTIONS] MATRIX';
  
  arguments = parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf));
  opt = arguments$options
  if (opt$verbose) print_options(opt)
  files = arguments$args
  out = opt$"output-file"
  n_reads = opt$'n-reads'
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }
  
  # input arguments
  fmatrix = files[1]
  
  # read files
  if (opt$verbose) write('Loading matrix...',stderr())
  x = as.matrix(read.table(fmatrix,row.names=1,check.names=F))
  if (is_full_matrix(x)==FALSE) x = MatrixInverseRotate45(x)
  if (sum(rownames(x)!=colnames(x))>0) { write(paste('Error: row names and column names should be identical!',sep=''),stderr()); quit(save='no') }
  features = as.matrix(read.table(opt$features,row.names=1,check.names=F))
  colnames(features) = c('effective-length','GC-content','mappability')

  # match rownames to feature loci
  i = match(rownames(x),rownames(features))
  n_unmatched = sum(is.na(i))
  if (n_unmatched>0) write(paste('Warning: ',n_unmatched,' loci were not matched to feature matrix!',sep=''),stderr())

  # set ignored loci rows/columns to zero
  ignored_rows = integer(0)
  ignored_cols = integer(0)
  if (opt$'ignored-loci'!="") {
    ignored_rows = sort(which(!is.na(match(rownames(x),as.vector(t(read.table(opt$'ignored-loci',check.names=F)))))))
    ignored_cols = sort(which(!is.na(match(colnames(x),as.vector(t(read.table(opt$'ignored-loci',check.names=F)))))))
    if (opt$verbose) { write(paste('Number of rows set to zero in input matrix (ignored loci) = ',length(ignored_rows),sep=''),stderr()); }
    if (opt$verbose) { write(paste('Number of columns set to zero in input matrix (ignored loci) = ',length(ignored_cols),sep=''),stderr()); }
    x[ignored_rows,] = 0          # ignored bins are set to zero!
    x[,ignored_cols] = 0
  }

  # normalize by effective length and total number of reads
  if (opt$"scale"==TRUE) {
    if (opt$verbose) write('Scaling matrix...',stderr())
    efflen = as.matrix(features[i,'effective-length'])
    efflen2 = efflen %*% t(efflen)
    mappability = as.matrix(features[i,'mappability'])
    mappability2 = mappability %*% t(mappability)
    y = x/(efflen2/1e6)/(n_reads/1e9)
    y[(efflen2<opt$'min-efflen'^2)|(mappability2<opt$'min-mappability'^2)] = NA         # regions with low effective length and/or low mappability will be considered as missing values
    if (opt$verbose) {
      write(paste('Max value of original matrix = ',max(x,na.rm=TRUE),sep=''),stderr())
      write(paste('Max value of scaled matrix = ',max(y,na.rm=TRUE),sep=''),stderr())
    }
  } else {
    y = x
  }

  # impute symmetric matrix
  if (opt$"impute"==TRUE) {
    if (opt$verbose) write('Imputing matrix...',stderr())
    y = MatrixImpute(y,ignored_rows,ignored_cols)
    if (opt$verbose) {
      write(paste('Symmetricity of original matrix = ',sum(x==t(x))/length(x),sep=''),stderr())
      write(paste('Symmetricity of imputed matrix = ',sum(y==t(y))/length(y),sep=''),stderr())
    }
  }
  
  if (opt$dist>0) {
    # rotate matrix
    if (opt$verbose) write('Rotating matrix...',stderr())
    y = MatrixRotate45(y,opt$dist)
  }
  
  # print new matrix
  if (opt$verbose) write("Saving new matrix...",stderr())
  write.table(format(y,scientific=TRUE,digits=4),row.names=TRUE,col.names=ifelse(opt$dist>0,FALSE,TRUE),quote=FALSE,sep='\t',file=out)

  # done
  if (opt$verbose) write("Done.",stderr())
  quit(save='no')
}


  



# ##########################
#  OPERATION = STANDARDIZE
# ##########################

op_standardize <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("--row-labels"), action="store_true",default=FALSE, help="Input matrix has row labels"),
    make_option(c("--stat-file"), default="", help="File containing Hi-C average counts as a function of distance [default=%default].")
  )
  usage = 'hic-matrix.r standardize [OPTIONS] MATRIX';
  
  arguments <- parse_args(args=cmdline_args,OptionParser(usage=usage,option_list=option_list),positional_arguments=c(0,Inf));
  opt <- arguments$options
  files <- arguments$args;
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }
  
  # input arguments
  fmatrix = files[1]
  row_labels = opt$"row-labels"
  fstat = opt$'stat-file'
  
  # read files
  if (opt$verbose) { write('Loading matrix...',stderr()); }
  if (row_labels==FALSE) { x = as.matrix(read.table(fmatrix,check.names=F)) } else { x = as.matrix(read.table(fmatrix,row.names=1,check.names=F)) }

  # standardize input matrix
  s = as.matrix(read.table(fstat,header=TRUE))[,1]
  max_dist = length(s)
  D = 1+abs(row(x)-col(x))
  y = x
  y[D>max_dist] = 0
  for (d in 1:max_dist) {
    z = y[D==d]
    y[D==d] = z/mean(z)*s[d]
  }

  # print new matrix
  if (opt$verbose) { write("Printing new matrix...",stderr()); }
  write.table(format(y,scientific=TRUE,digits=4),quote=F,sep='\t')

  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no')
}


  



# ########################
#  OPERATION = SUBMATRIX
# ########################

op_submatrix <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("--start"), default=1, help="Start row/column [default \"%default\"]."),
    make_option(c("--end"), default=1, help="End row/column [default \"%default\"].")
  );
  usage = 'hic-matrix.r submatrix [OPTIONS] MATRIX';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  files <- arguments$args;
  a = as.integer(opt$start)
  b = as.integer(opt$end)
  
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input file
  f <- files[1]

  # load matrices  
  if (opt$verbose) { write("Loading matrix...",stderr()); }
  m = as.matrix(read.table(f,check.names=F))

  # create submatrix
  if (opt$verbose) { write("Printing submatrix...",stderr()); }
  m = m[a:b,a:b]
  write.table(m,quote=F,sep='\t')

  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}


  



# ########################
#  OPERATION = TRANSLOC
# ########################

op_transloc <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("--dist"), default=0, help="Distance in number of bins [default \"%default\"].")
  );
  usage = 'hic-matrix.r transloc [OPTIONS] MATRIX';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  files <- arguments$args;
  dist = as.integer(opt$dist)
  
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input file
  f <- files[1]

  # load matrices  
  if (opt$verbose) { write("Loading matrix...",stderr()); }
  m = as.matrix(read.table(f,row.names=1,header=FALSE,check.names=F))

  # compute translocation scores
  if (opt$verbose) { write("Computing translocation scores...",stderr()); }
  m = log2(m+1)
  f = function(x,y) { mean(m[x:(x+dist-1),(y-dist+1):y]) }
  v = rep(0,(nrow(m)-dist)^2)
  k = 1
  for (a in 7000:(7300-dist)) { 
    print(a)
    for (b in (14500+dist):(nrow(m)-dist)) { 
      p = f(a,a)+f(b,b)
      q = f(a,b)+f(b,a)
      v[k] = q/(p+q) 
      k = k + 1
    } 
  }
  v = v[1:(k-1)]
  
#  write.table(m,quote=F,sep='\t')

  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}


  



# ########################
#  OPERATION = STATS
# ########################

op_stats <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
    make_option(c("-b","--bin-size"), default=0, help="Bin size in nucleotides [default \"%default\"]."),
#    make_option(c("-L","--sample-labels"), default="", help="Comma-separated list of sample labels (required) [default \"%default\"]."),
    make_option(c("--features"), default="", help="Feature file [default \"%default\"].")
  );
  usage = 'hic-matrix.r stats [OPTIONS] MATRIX-DIRECTORIES';
  
  # get command line options (if help option encountered print help and exit)
  arguments = parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
  opt = arguments$options
  inp_dirs = arguments$args
  out_dir = opt$'output-dir'
  bin_size = as.integer(opt$"bin-size")
#  if (opt$"sample-labels"=="") { write('Error: please specify sample labels!',stderr()); quit(save='no'); }
#  sample_labels = unlist(strsplit(opt$"sample-labels",split=','))
  sample_labels = sub(".*/","",inp_dirs)

  # create output directory
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(out_dir)==FALSE) dir.create(out_dir)

  # load matrices
  n_bins = 50000000/bin_size
  D = matrix(0,length(inp_dirs),n_bins)
  k = 1
  while (k <= length(inp_dirs)) {
    files = Sys.glob(paste(inp_dirs[k],"matrix.*.tsv",sep="/"))
    for (f in files) {
      if (opt$verbose) write(paste("Processing matrix '",f,"'...",sep=''),stderr())
      x = as.matrix(read.table(f,check.names=F))
      x[is.na(x)] = 0
      counts = calc_counts_per_distance(x)                     # average Hi-C count per distance for each matrix
      counts = c(counts,rep(0,n_bins))[1:n_bins]               # truncate or pad with zeroes
      D[k,] = D[k,] + counts
    }
    D[k,] = D[k,]/length(files)
    k = k + 1
  }

  # create Hi-C count vs distance plot
  pdf(paste(out_dir,'/stats.pdf',sep=''))
  colors = c('blue','green4','red','magenta','yellow','cyan','brown','pink','orange','gray','black','purple')
  colors = rep(colors,ceiling(nrow(D)/length(colors)))[1:nrow(D)]
  par(mfrow=c(2,1),mar=c(2,4,2,2))
  ylim = c(min(D[D>0]),max(D))
  d = (1:ncol(D))*bin_size/1000
  plot(d,D[1,],type='l',col=colors[1],log="xy",xlab='distance (kb)',ylab='Average normalized Hi-C count',main='Hi-C count as a function of distance',ylim=ylim)
  for (k in 1:nrow(D)) lines(d,D[k,],col=colors[k])
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x='top',cex=0.70,legend=sample_labels,col=colors,lwd=2)
  dev.off()  

  quit(save='no')

  # correlations with features
  if (opt$features!="") features = as.matrix(read.table(opt$features,row.names=1,check.names=F))
  if (opt$features!="") {
    if (opt$verbose) { write("Computing correlations with features...",stderr()); }
    a = which(rownames(features)==rownames(y)[1])
    b = which(rownames(features)==rownames(y)[nrow(y)])
    features = features[a:b,]
    I = abs(row(y)-col(y))
    for (k in 1:ncol(features)) {
      par(mfrow=c(2,2))
      e = features[,k] %*% t(features[,k])
      dist = as.integer(c(0.1,0.2,0.5)*nrow(y))
      dist = unique(c(1,dist))
      for (d in dist) {
        J = I<d 
        ee = as.vector(e[J])
        x = sapply(1:nrow(solObj),function(k) cor(ee,as.vector(solObj[k,,][J])))
        plot(x,type='l',col='red',xlab='parameter',ylab='correlation',main=paste('d=',d,sep=''))
      }
    }
  }
  
  # done
  save(files,D,file=paste(out_dir,'/stats.RData',sep=''))
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no')

}


  




# ########################
#  OPERATION = DIFF
# ########################

op_loopdiff <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("--debug"), action="store_true",default=FALSE, help="Debugging mode (using matrix inverse-rotate)."),
    make_option(c("--lambda-id"), default=1, help="Lambda to be selected from RData file (if applicable) [default \"%default\"]."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
    make_option(c("-b","--bin-size"), default=0, help="Bin size in nucleotides [default \"%default\"]."),
    make_option(c("--n-reads"), default="", help="Comma-separated list of total number of reads (required) [default \"%default\"]."),
    make_option(c("--min-count"), default=20, help="Minimum Hi-C count (scaled) [default \"%default\"]."),
    make_option(c("--loop-cutoff"), default=5.0, help="Loop enrichment score cutoff [default \"%default\"]."),
    make_option(c("--min-distance"), default=100000, help="Minimum loop distance in nucleotides [default \"%default\"]."),
    make_option(c("--bins"), default="", help="Comma-separated list of bins to be included in the snapshot [default \"%default\"]."),
    make_option(c("-L","--sample-labels"), default="", help="Comma-separated list of sample labels (required) [default \"%default\"].")
  );
  usage = 'hic-matrix.r diff [OPTIONS] MATRIX-1(tsv/RData) MATRIX-2(tsv/RData)';
  
  # get command line options (if help option encountered print help and exit)
  arguments = parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt = arguments$options
  files = arguments$args
  if (length(files)!=2) { write('Error: this operation requires two input matrices!',stderr()); quit(save='no'); }

  # input arguments
  out_dir = opt$'output-dir'
  n_reads = as.integer(strsplit(opt$'n-reads',',')[[1]])
  bin_size = opt$'bin-size'
  min_count = opt$'min-count'
  loop_cutoff = opt$'loop-cutoff'
  min_dist = opt$'min-distance'
  sample_labels = opt$'sample-labels'
  bins = as.integer(strsplit(opt$bins,',')[[1]])
  if (sample_labels=="") { write('Error: please specify sample labels!',stderr()); quit(save='no'); }
  sample_labels = strsplit(sample_labels,split=',')[[1]]

  # Create output directory
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Error: output directory already exists!',stderr()); quit(save='no'); }
  if (bin_size<=0) { write('Error: please specify a valid bin size!',stderr()); quit(save='no'); }
  if (length(n_reads)==0) n_reads = c(1,1)       # no scaling
  if (length(n_reads)!=2) { write('Error: please provide the total number of reads for all samples!',stderr()); quit(save='no'); }

  # load matrices  
  if (opt$verbose) { write("Loading matrices...",stderr()); }
  m = { }
  m_scaled = { }
  m_norm = { }
  m_zscore = { }
  for (f in 1:length(files)) {
    ext = strsplit(files[f],'.*[.]')[[1]][2]
    if (ext=='RData') {
      env1 = new.env()
      load(files[f],env1)
      if ((env1$opt$preprocess!='dist')&&(env1$opt$preprocess!='distlog2')) { write("Error: estimation preprocessing should be dist or distlog2 for this operation!",stderr()); quit(save='no') }
      if ((opt$'lambda-id'<1)||(opt$'lambda-id'>length(env1$lambdas))) { write("Error: lambda id out of bounds!",stderr()); quit(save='no') }
      m[[f]] = env1$x                                                                                   # raw matrix
      m_norm[[f]] = env1$solObj[opt$'lambda-id',,]                                                      # estimated matrix
    } else {
      m[[f]] = as.matrix(read.table(files[f],check.names=F,row.names=1))                                # row labels are necessary: should include locus information
      if ((opt$debug==TRUE)&(is_full_matrix(m[[f]])==FALSE)) m[[f]] = MatrixInverseRotate45(m[[f]])     # only for debugging
      m_norm[[f]] = PreprocessMatrix(m[[f]],'dist',pseudo=1,cutoff=0)
    }
    m_scaled[[f]] = m[[f]]/(n_reads[f]/max(n_reads))                                                    # scaled matrix
    m_zscore[[f]] = PreprocessMatrix(m[[f]],'zscore',pseudo=1,cutoff=0)                                 # z-scores
  }
  if (sum(apply(sapply(m,dim),1,function(v) {max(v)-min(v)})>0)!=0) { write('Error: input matrices should have the same dimensions!',stderr()); quit(save='no'); }
  full_matrix = is_full_matrix(m[[1]])

  # Filtering: (a) low Hi-C scores, (b) low observed/expected or (c) small distance
  if (full_matrix) { D = col(m[[1]])-row(m[[1]])
  } else { D = col(m[[1]])-1 }
  D = D*bin_size/1000
  i_removed = ((m_scaled[[1]]<min_count)&(m_scaled[[2]]<min_count)) | ((m_norm[[1]]<loop_cutoff)&(m_norm[[2]]<loop_cutoff)) | (abs(D)<min_dist/1000)
  if (sum(!i_removed)==0) { write("No interactions found!", stderr()); quit(save='no') }
  
  # create plots
  pdf(paste(out_dir,'/diff.pdf',sep=''))

  # calculate log2 fold-changes
  if (opt$verbose) { write("Calculating log2 fold-changes...",stderr()); }
  pseudo_scaled = min_count/2
  pseudo_norm = 0.5
  fold_scaled = log2((m_scaled[[2]]+pseudo_scaled)/(m_scaled[[1]]+pseudo_scaled))
  fold_norm = log2((m_norm[[2]]+pseudo_norm)/(m_norm[[1]]+pseudo_norm))
  diff_zscore = m_zscore[[2]]-m_zscore[[1]]

  # plot heatmaps and histogram
  if (opt$verbose) write("Creating plots...",stderr())
  par(mfrow=c(2,2),cex=0.75)

  # scores normalized by total reads
  x = as.vector(m_scaled[[1]][!i_removed])
  y = as.vector(m_scaled[[2]][!i_removed])
  smoothScatter(log2(x),log2(y),xlab=sample_labels[1],ylab=sample_labels[2],main='Hi-C scaled counts')
  hist(fold_scaled[!i_removed],main='Histogram of log2 fold-changes',xlab='log2 fold-change',ylab='Frequency',n=20)

  # distance-normalized scores
  x = as.vector(m_norm[[1]][!i_removed])
  y = as.vector(m_norm[[2]][!i_removed])
  smoothScatter(log2(x),log2(y),xlab=sample_labels[1],ylab=sample_labels[2],main='Hi-C distance-normalized scores')
  hist(fold_norm[!i_removed],main='Histogram of log2 fold-changes',xlab='log2 fold-change',ylab='Frequency',n=20)

  # z-scores
  x = as.vector(m_zscore[[1]][!i_removed])
  y = as.vector(m_zscore[[2]][!i_removed])
  smoothScatter(x,y,xlab=sample_labels[1],ylab=sample_labels[2],main='Hi-C distance-normalized z-scores')
  hist(diff_zscore[!i_removed],main='Histogram of log2 fold-changes',xlab='log2 fold-change',ylab='Frequency',n=20)

  # heatmaps
  if (full_matrix) {
    n = nrow(m[[1]])
    cc = 0
    xpos = fold_scaled; xpos[xpos<=cc] = NA
    xneg = -fold_scaled; xneg[xneg<=cc] = NA
    if (length(bins)>0) {
      i_sampled = max(1,min(bins)):min(n,max(bins))
      bin_points = (bins-min(bins))/(max(bins)-min(bins))
      bin_points = bin_points[-length(bin_points)]
      bin_points = bin_points[-1]
    } else {
      i_sampled = seq(1,n,length.out=200)
      bin_points = c()
    }
#    image(max(xpos,na.rm=TRUE)-xpos[i_sampled,i_sampled],col='red',main='Hi-C gain')
    image(xpos,col='red',main='Hi-C gain')
    for (pp in bin_points) draw.circle(pp,pp,0.01,col='green4',border='green4')
#    image(max(xneg,na.rm=TRUE)-xneg[i_sampled,i_sampled],col='blue',main='Hi-C loss')
    image(xneg,col='blue',main='Hi-C loss')
    for (pp in bin_points) draw.circle(pp,pp,0.01,col='green4',border='green4')
  }
  
  # plot is done
  dev.off()

  # Save fold changes in matrix format
#  if (opt$verbose) { write("Storing results...",stderr()); }
#  write.table(round(fold_scaled,4),quote=F,row.names=TRUE,col.names=full_matrix,sep='\t',file=paste(out_dir,'/matrix_scaled.tsv',sep=''))
#  write.table(round(fold_norm,4),quote=F,row.names=TRUE,col.names=full_matrix,sep='\t',file=paste(out_dir,'/matrix_norm.tsv',sep=''))
#  write.table(round(diff_zscore,4),quote=F,row.names=TRUE,col.names=full_matrix,sep='\t',file=paste(out_dir,'/matrix_zscore.tsv',sep=''))

  # Save fold changes in table format
  if (full_matrix) {
    labels = melt(fold_scaled)[!i_removed,1:2]
  } else {
    rnames = rownames(m[[1]])
    N = length(rnames)
    f <- function(pos) { i = (pos-1)%%N+1; j = (pos-1)%/%N+1; x = i-(j-1)%/%2; y = x+j-1; return(c(rnames[x],rnames[y])) }
    labels = t(sapply(which(!i_removed),f))
  }
  table = cbind(labels,round(fold_scaled[!i_removed],3))
  table = cbind(table,round(fold_norm[!i_removed],3))
  table = cbind(table,round(diff_zscore[!i_removed],3))
  for (z in c(m,m_scaled,m_norm,m_zscore)) table = cbind(table,round(z[!i_removed],3))
  table = cbind(table,as.integer(1000*D[!i_removed]))
  colnames(table) = c('locus1','locus2','log2-fold-scaled','log2-fold-norm','diff-zscore','sample1-raw','sample2-raw','sample1-scaled','sample2-scaled','sample1-norm','sample2-norm','sample1-zscore','sample2-zscore','distance')
  # add symmetric interactions for distance-restricted matrix
  if (full_matrix==FALSE) {  
    table2 = table[,c(2,1,3:ncol(table))]
    colnames(table2) = colnames(table)
    table2[,'distance'] = -as.integer(table2[,'distance'])
    table = rbind(table,table2)
  }
  write.table(table,file=paste(out_dir,'/table.tsv',sep=''),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
  
  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no')
}







# ########################
#  OPERATION = LOOPS
# ########################

op_loops <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
    make_option(c("-L","--sample-labels"), default="", help="Comma-separated list of sample labels (required) [default \"%default\"]."),
    make_option(c("--n-reads"), default="", help="Comma-separated list of total number of reads per sample (required) [default \"%default\"]."),
    make_option(c("--bin-size"), default=0, help="Bin size in nucleotides for RPKB calculation [default \"%default\"]."),
    make_option(c("--lambda-id"), default=1, help="Lambda to be selected from RData file (if applicable) [default \"%default\"]."),
    make_option(c("--rpk2b-cutoff"), default=1.0, help="RPK2B score cutoff [default \"%default\"]."),
    make_option(c("--loop-cutoff"), default=5.0, help="Distance-normalized loop score cutoff [default \"%default\"]."),
    make_option(c("--min-distance"), default=40000, help="Minimum loop distance in nucleotides [default \"%default\"]."),
    make_option(c("--debug"), action="store_true",default=FALSE, help="Debugging mode (using matrix inverse-rotate).")
  )
  usage = '\
\thic-matrix.r loops [OPTIONS] MATRIX(tsv/RData) ...\
\
Function: \
\tIdentify Hi-C interactions in input contact matrices. \
\
Input: \
\tMultiple compatible matrices in tsv or RData format'
  
  # Get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
  opt <- arguments$options
  if (opt$verbose) print_options(opt)
  files <- arguments$args
  if (length(files)==0) { write(paste('Usage:',usage),stderr()); quit(save='no') }

  # Input arguments
  out_dir = opt$'output-dir'
  n_reads = as.integer(strsplit(opt$'n-reads', ",")[[1]])
  bin_size = opt$'bin-size'
  rpk2b_cutoff = opt$'rpk2b-cutoff'
  loop_cutoff = opt$'loop-cutoff'
  min_dist = opt$'min-distance'
  sample_labels = opt$'sample-labels'

  # Error checking
  n_files = length(files)
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Warning: output directory already exists, overwriting...',stderr()) }
  if (length(n_reads)==0) n_reads = rep(0,n_files)       # no scaling
  if (length(n_reads)!=n_files) { write('Error: please provide the total number of reads for all samples!',stderr()); quit(save='no'); }
  if (bin_size<=0) { write('Error: please specify a valid bin size!',stderr()); quit(save='no'); }
  if (sample_labels=="") { write('Error: please specify sample labels!',stderr()); quit(save='no'); }
  sample_labels = strsplit(sample_labels,split=',')[[1]]
  if (length(sample_labels)!=n_files) { write('Error: please specify sample labels for each sample!',stderr()); quit(save='no'); }

  # Load all matrices
  full_matrix = rep(FALSE,n_files)
  x_raw = {}
  x_scaled = {}
  x_norm = {}
  prep = 'dist'  #'zscore' # 'dist'
  for (f in 1:n_files) {
    if (opt$verbose) write(paste("Loading input matrix '",files[f],"'...",sep=''),stderr())
    ext = strsplit(files[f],'.*[.]')[[1]][2]
    if (ext=='RData') {
      est = new.env()
      load(files[f],est)
      #if ((est$opt$preprocess!='zscore')&&(est$opt$preprocess!='dist')&&(est$opt$preprocess!='distlog2')) { write("Error: estimation preprocessing should be dist or distlog2 for this operation!",stderr()); quit(save='no') }
      if ((opt$'lambda-id'<1)||(opt$'lambda-id'>length(est$lambdas))) { write("Error: lambda id out of bounds!",stderr()); quit(save='no') }
      x_raw[[f]] = est$x                                                                                            # raw matrix
      x_norm[[f]] = est$solObj[opt$'lambda-id',,]                                                                   # estimated distance-normalized matrix
    } else {  
      x_raw[[f]] = as.matrix(read.table(files[f],check.names=F,row.names=1))                                        # row labels are necessary (locus information)
      if ((opt$debug==TRUE)&(is_full_matrix(x_raw[[f]])==FALSE)) x_raw[[f]] = MatrixInverseRotate45(x_raw[[f]])     # only for debugging
      x_norm[[f]] = PreprocessMatrix(x_raw[[f]],prep,pseudo=1,cutoff=0)                                             # distance-normalized matrix
    }
    a = ifelse(n_reads[f]<=0, 1, (n_reads[f]/1e9)*(bin_size/1000)^2)
    x_scaled[[f]] = x_raw[[f]]/a                                                                                    # scaled matrix
    full_matrix[f] = is_full_matrix(x_raw[[f]])
  }
  
  # check matrix compatibility
  if ((length(unique(full_matrix))>1)||(sum(apply(sapply(x_raw,dim),1,function(v) length(unique(v))!=1))>0)) { write("Error: input matrices are not compatible!\n",stderr()); quit(save='no') }

  # Filtering: (a) low Hi-C scores, (b) low observed/expected or (c) small distance
  if (opt$verbose) write("Filtering interactions...",stderr())
  if (full_matrix[1]) { D = col(x_raw[[1]])-row(x_raw[[1]]) } else { D = col(x_raw[[1]])-1 }
  D = D*bin_size
  z = sapply(1:n_files,function(k) (x_norm[[k]]>=loop_cutoff) & (x_scaled[[k]]>=rpk2b_cutoff))         # apply scaled count and distance-normalized score cutoffs
  i_selected = apply(z,1,max) & (abs(as.vector(D))>=min_dist)                                          # apply minimum distance cutoff
  n_selected = sum(i_selected)
  if (opt$verbose) write(paste("Selected interactions = ",n_selected,sep=''),stderr())
  if (n_selected==0) { write("No interactions found!", stderr()); system(paste('touch ',out_dir,'/loops.tsv',sep='')); quit(save='no') }
  
  # Generate table
  if (opt$verbose) write("Generating table...",stderr())
  if (full_matrix[1]==TRUE) {                                    # get loci labels
    table = cbind(melt(x_raw[[1]]))[i_selected,1:2,drop=FALSE]
  } else {
    rnames = rownames(x_raw[[1]])
    N = length(rnames)
    get_labels <- function(pos) { i = (pos-1)%%N+1; j = (pos-1)%/%N+1; x = i-(j-1)%/%2; y = x+j-1; return(c(rnames[x],rnames[y])) }
    table = t(sapply(which(i_selected),get_labels))
  }
  for (v in c(x_raw,x_scaled,x_norm)) table = cbind(table,round(as.vector(v)[i_selected],3))            # merge data into a single table
  table = cbind(table,as.integer(D[i_selected]))
  colnames(table) = c('locus1','locus2',sapply(sample_labels,paste,'-hic-count',sep=''),sapply(sample_labels,paste,'-hic-scaled',sep=''),sapply(sample_labels,paste,'-hic-distnorm',sep=''),'distance')
  if (full_matrix[1]==FALSE) {                                   # reverse order of loci to account for loop symmetricity (not necessary for full matrices)
    table2 = table[,c(2,1,3:ncol(table)),drop=FALSE]
    colnames(table2) = colnames(table)
    table2[,'distance'] = -as.integer(table2[,'distance'])
    table = rbind(table,table2)
  }

  # Save results
  if (opt$verbose) write("Saving results...",stderr())
  write.table(table,file=paste(out_dir,'/loops.tsv',sep=''),row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
  
  # Create plots
  if (opt$verbose) write("Creating plots...",stderr())
  pdf(paste(out_dir,'/plots.pdf',sep=''))
  par(mfrow=c(2,2))
  boxplot(lapply(x_scaled,function(mat) mat[mat>=rpk2b_cutoff]),log='y',names=sample_labels,main='Scaled Hi-C counts')
  boxplot(lapply(x_norm,function(mat) mat[mat>=loop_cutoff]),log='y',names=sample_labels,main='Distance-normalized Hi-C scores')
  counts = lapply(x_scaled,function(mat) { z=calc_counts_per_distance(mat); z=z/max(z); return(z)})
  colors = c('red','green','blue','magenta','black','purple','yellow')
  colors = as.vector(replicate(length(counts)%/%length(colors)+1,colors))[1:length(counts)]
  x = 1:length(counts[[1]])*bin_size/1000
  ylim = c(min(sapply(counts,function(v) min(v[v>0]))),max(sapply(counts,max)))
  plot(x,counts[[1]],type='l',col='red',log="xy",xlab='distance (kb)',ylab='Scaled Hi-C count (RPK2B)',main='Hi-C scaled count per distance',ylim=ylim)
  for (k in 1:length(counts)) lines(x,counts[[k]],col=colors[k])
  legend('bottomleft',sample_labels,cex=0.75,pch=16,col=colors,inset=0.05);
#  plot(abs(D[i]),y[i],log='xy',pch=18,col='green4',xlab='distance in kb (log scale)',ylab='Loop enrichment score (log scale)',main='Loop score as a function of distance')
#  hist(abs(D[i&(abs(D)<10000000)]),n=100,xlab='distance in kb',main='Histogram of loops per distance')
  dev.off()
    
  # done
  if (opt$verbose) write("Done.",stderr())
  quit(save='no')
}


  





# ########################
#  OPERATION = MATRICES
# ########################

op_matrices <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
    make_option(c("--min-lambda"), default=0.0, help="Minimum value for parameter lambda [default=%default]."),
    make_option(c("--max-lambda"), default=1.0, help="Maximum value for parameter lambda [default=%default]."),
    make_option(c("--n-lambda"), default=6, help="Number of lambdas [default=%default]."),
    make_option(c("--log2-lambda"), action="store_true",default=FALSE, help="Use log2 scale for lambda range."),
    make_option(c("--gamma"), default=0.0, help="Value for sparsity parameter gamma [default=%default].")
  )
  usage = 'hic-matrix.r matrices [OPTIONS] ESTIMATED-RDATA-FILE';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  files <- arguments$args;
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  f <- files[1];
  out_dir <- opt$'output-dir';

  # create output directory
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Warning: output directory already exists!',stderr()) }

  # load estimated matrices (RData file)
  if (opt$verbose) write("Loading data...",stderr())
  e = new.env()
  load(f,e)
  
  # inverse rotate45 (TODO: make this an option?)
  invrotate = is_full_matrix(e$x)==FALSE
  
  # save matrices in text format
  if (e$opt$"max-lambda"==Inf) {
    # create lambda values (zero always included if log2 scale)
    if (opt$"log2-lambda"==FALSE) { lambdas = unique(seq(opt$"min-lambda",opt$"max-lambda",length.out=opt$"n-lambda"))
    } else { lambdas = unique(c(0,2^seq(log2(max(opt$"min-lambda",0.01)),log2(opt$"max-lambda"),length.out=opt$"n-lambda"-1))) }

    # create new estimated matrices for specific lambdas/gammas
    if (opt$verbose) write('Saving estimated matrices for specific lambdas/gammas...',stderr())
    n_matrices = length(lambdas)
    for (k in 1:n_matrices) {
      z = GetSolution(e$solObj,n_rows=nrow(e$y),invrotate=invrotate,lambda=lambdas[k],gamma=opt$gamma)
      rownames(z) = rownames(e$y)
      colnames(z) = colnames(e$y)
      fout = paste(out_dir,'/matrix.k=',formatC(k,width=3,format='d',flag='0'),'.tsv',sep='')
      write.table(format(z,scientific=TRUE,digits=4),row.names=TRUE,col.names=is_full_matrix(z),quote=FALSE,sep='\t',file=fout)
    }
  
  } else {
    # save matrices for each lambda value
    if (opt$verbose) write("Saving estimated matrices...",stderr())
    n_matrices = dim(e$solObj)[1]
    for (k in 1:n_matrices) {
      z = e$solObj[k,,]
      if (invrotate==TRUE) z = MatrixInverseRotate45(z)
      rownames(z) = rownames(e$y)
      colnames(z) = colnames(e$y)
      fout = paste(out_dir,'/matrix.k=',formatC(k,width=3,format='d',flag='0'),'.tsv',sep='')
      write.table(format(z,scientific=TRUE,digits=4),row.names=TRUE,col.names=is_full_matrix(z),quote=FALSE,sep='\t',file=fout)
    }
  }

  # done
  if (opt$verbose) write("Done.",stderr())
  quit(save='no')
}


  




# ########################
#  OPERATION = HEATMAPS
# ########################

op_heatmaps <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-d","--heatmap-size"), default=600, help="Heatmap image size [default \"%default\"]."),
    make_option(c("-o","--output-file"), default="", help="Output pdf file (required) [default \"%default\"].")
  );
  usage = 'hic-matrix.r heatmaps [OPTIONS] ESTIMATED-RDATA-FILE';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  files <- arguments$args;
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  f <- files[1];
  heatmap_size <- opt$"heatmap-size";
  out <- opt$"output-file";

  # create matrices  
  if (opt$verbose) { write("Loading data...",stderr()); }
  load(f);
  if (opt$verbose) { write("Plotting contact maps...",stderr()); }
  pdf(out)
  n_matrices = dim(solObj)[1]
  n = dim(solObj)[2]
  I = as.integer(seq(1,n,length.out=min(n,heatmap_size)));   # sample matrix for heatmap generation
  par(mfrow=c(2,2))
  for (k in 1:n_matrices) {
    x = solObj[k,I,I]
    xpos = x; xpos[xpos<0] = 0; 
    xneg = -x; xneg[xneg<0] = 0;
    image(xpos,main=paste('lambda=',lambdas[k],sep=''))
    image(xneg,main=paste('lambda=',lambdas[k],sep=''))
  }
  dev.off()
  
  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}


  




# ########################
#  OPERATION = SNAPSHOTS
# ########################

op_snapshots <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-u","--upper-bound"), default=Inf, help="Upper bound on matrix values [default \"%default\"]."),
    make_option(c("-l","--lower-bound"), default=0, help="Lower bound on matrix values [default \"%default\"]."),
    make_option(c("-d","--max-distance"), default=100, help="Max off-diagonal distance in bins [default \"%default\"]."),
    make_option(c("-b","--bins"), default="1,1000", help="Comma-separated list of bins to be included in the snapshot [default \"%default\"]."),
    make_option(c("-o","--output-file"), default="", help="Output pdf file (required) [default \"%default\"].")
  );
  usage = 'hic-matrix.r snapshots [OPTIONS] ESTIMATED-RDATA-FILE(S)';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  files <- arguments$args;
  if (length(files)<1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  max_dist = as.integer(opt$"max-distance")
  bins = as.integer(strsplit(opt$bins,',')[[1]])
  out = opt$"output-file"
  U = as.double(opt$'upper-bound')
  L = as.double(opt$'lower-bound')
  
  # create matrices  
  pdf(out)
  n_panels = 7
  for (f in 1:length(files)) { 
    if (opt$verbose) { write(paste("Loading file ",files[f],"...",sep=''),stderr()); }
    load(files[f])
    if (opt$verbose) { write("Creating snapshots...",stderr()); }
    full_matrix = is_full_matrix(x)
    n_matrices = dim(solObj)[1]
    n = dim(solObj)[2]
    I = max(1,min(bins)):min(n,max(bins))
    bb = (bins-min(I))/(max(I)-min(I))
    bb = bb[(bb>0)&(bb<1)] 
    par(mfrow=c(n_panels,1),mar=c(1,1,1,1))
    if (full_matrix) { z = MatrixRotate45(t(x[I,I]),max_dist) } else { z = x[I,] }        # raw contact matrix
    z = log2(z+1)
    z[z<=L] = L; z[z>=U] = U; if (U<=0) z=-z; z = (z-min(z))/(max(z)-min(z))
    image(1-z,main='raw matrix',xaxt='n',yaxt='n')
    for (j in 1:length(bb)) draw.circle(bb[j],0,0.005,col='blue',border='blue')
    if (full_matrix) { z = MatrixRotate45(t(y[I,I]),max_dist) } else { z = y[I,] }         # normalize matrix (gamma=0, lambda=0)
    z[z<=L] = L; z[z>=U] = U; if (U<=0) z=-z; z = (z-min(z))/(max(z)-min(z))
    image(1-z,main=substitute(paste(gamma,'=0, ',lambda,'=0',sep='')),xaxt='n',yaxt='n')
    for (j in 1:length(bb)) draw.circle(bb[j],0,0.005,col='blue',border='blue')
    for (k in seq(1,n_matrices,length.out=n_panels-2)) {
      if (full_matrix) { z = MatrixRotate45(t(solObj[k,I,I]),max_dist) } else { z = solObj[k,I,] }
      z[z<=L] = L; z[z>=U] = U; if (U<=0) z=-z; z = (z-min(z))/(max(z)-min(z))
      image(1-z,main=substitute(paste(gamma,'=',GAMMA,', ',lambda,'=',LAMBDA,sep=''),list(GAMMA=gammas[1],LAMBDA=round(lambdas[k],3))),xaxt='n',yaxt='n')
      for (j in 1:length(bb)) draw.circle(bb[j],0,0.005,col='blue',border='blue')
    }
  }
  dev.off()
  
  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}


  

###### ExtractEstimations

ExtractEstimations = function(filename, opt)
{ 
  # load estimation
  e = new.env()
  load(filename,e)

  if (e$opt$"max-lambda"==Inf) {
    # setup
    algorithm = e$opt$algorithm
    n = nrow(e$x)
    m = ncol(e$x)
  
    # create lambdas
    e$opt$"max-lambda" = opt$"max-lambda"
    e$opt$"min-lambda" = opt$"min-lambda"
    e$opt$"n-lambda" = opt$"n-lambda"
    e$opt$"log2-lambda" = opt$"log2-lambda"
    e$opt$gamma = opt$gamma
    if (opt$"max-lambda"==Inf) { lambdas = NULL;   # TODO: Error message!
    } else if (opt$"log2-lambda"==FALSE) { lambdas = unique(seq(opt$"min-lambda",opt$"max-lambda",length.out=opt$"n-lambda"));
    } else { lambdas = unique(c(0,2^seq(log2(max(opt$"min-lambda",0.01)),log2(opt$"max-lambda"),length.out=opt$"n-lambda"-1))); }
    n_matrices = length(lambdas)
    e$lambdas = lambdas
    e$gammas = opt$gamma
    
    # create new estimated matrices for specific lambdas/gammas
    if (opt$verbose) { write('Extracting solutions for specific lambdas/gammas...',stderr()); }
    solObj0 = e$solObj
    e$solObj = array(0,dim=c(n_matrices,n,m))
    for (i in 1:n_matrices) e$solObj[i,,] = GetSolution(solObj0,n_rows=n,invrotate=FALSE,lambda=lambdas[i],gamma=opt$gamma)

  } else {
    write("Warning: specific lambdas/gammas are ignored (only applicable if estimated max-lambda was set to Inf)!",stderr())    

  }
  
  return(e)
}



###### LoadEstimation

LoadEstimation = function(filename, options, replace.na)
{
  est = {}
  ext = strsplit(filename,'.*[.]')[[1]][2]
  if (ext=="RData") {                                 # get estimated matrices from RData file
    load(filename)  
    est$x = x
    est$y = y
    est$ignored_rows = sort(unique(ignored_rows))
    est$ignored_cols = sort(unique(ignored_cols))     # TODO: ignored_cols is set to NULL, produces warning
    if (opt$"max-lambda"==Inf) {                      # extract specific lambdas (found in options) if max-lambda=Inf
      # create gamma/lambda values (zero always included if log2 scale)
      est$gammas = opt$gammas
      if (options$"log2-lambda"==FALSE) { lambdas = unique(seq(options$"min-lambda",options$"max-lambda",length.out=options$"n-lambda"))
      } else { lambdas = unique(c(0,2^seq(log2(max(options$"min-lambda",0.01)),log2(options$"max-lambda"),length.out=options$"n-lambda"-1))) }
      est$lambdas = lambdas
      
      # create new estimated matrices for specific lambdas/gammas
      if (options$verbose) write('Loading matrices estimated for specific lambdas/gammas...',stderr())
      n_matrices = length(lambdas)
      est$solObj = array(0,dim=c(n_matrices,dim(y)))
      for (k in 1:n_matrices) est$solObj[k,,] = GetSolution(solObj,n_rows=nrow(y),invrotate=FALSE,lambda=lambdas[k],gamma=options$gamma)

    } else {                                         # keep original lambdas
      est$gammas = gammas
      est$lambdas = lambdas
      est$solObj = solObj
    }
    
  } else {                                            # if input is tsv, then do all necessary formatting and preprocessing
    # read matrix from text file
    if (options$'row-labels'==TRUE) { est$x = as.matrix(read.table(filename,check.names=F,row.names=1))
    } else { est$x = as.matrix(read.table(filename,check.names=F,row.names=NULL)) }

    # set lambdas/gammas
    est$gammas = 0
    est$lambdas = 0
        
    # replace NAs with 0
    if (replace.na==TRUE) est$x[is.na(est$x)] = 0

    # set ignored rows/columns to zero
    est$ignored_rows = integer(0)
    est$ignored_cols = integer(0)
    fignored = options$"ignored-loci"
    if (fignored!="") {
      est$ignored_rows = sort(which(!is.na(match(rownames(est$x),as.vector(t(read.table(fignored,check.names=F)))))))
      est$ignored_cols = sort(which(!is.na(match(colnames(est$x),as.vector(t(read.table(fignored,check.names=F)))))))
      if (options$verbose) { write(paste('Number of rows set to zero in input matrix (ignored loci) = ',length(est$ignored_rows),sep=''),stderr()); }
      if (options$verbose) { write(paste('Number of columns set to zero in input matrix (ignored loci) = ',length(est$ignored_cols),sep=''),stderr()); }
      est$x[est$ignored_rows,] = 0
      est$x[,est$ignored_cols] = 0
    }

    # preprocess input matrix
    est$y = PreprocessMatrix(est$x,preprocess=options$preprocess,pseudo=1,cutoff=1)
    est$solObj = array(0,dim=c(1,dim(est$x)[1],dim(est$x)[2]))
    est$solObj[1,,] = est$y
  }

  # create submatrix
  if (options$bins!='') {
    bins = as.integer(strsplit(options$'bins',',')[[1]])
    n_matrices = dim(est$solObj)[1]
    n = nrow(est$y)
    bin_range = max(1,min(bins)):min(n,max(bins))
    est$ignored_rows = integer(0)
    est$ignored_cols = integer(0)
    solObj0 = est$solObj
    est$x = est$x[bin_range,bin_range]        # TODO: adapt for distance-restricted matrices
    est$y = est$y[bin_range,bin_range]
    n = length(bin_range)
    est$solObj = array(0,dim=c(n_matrices,n,n))
    for (k in 1:n_matrices) est$solObj[k,,] = solObj0[k,bin_range,bin_range]
  }

  return(est)
}




###### IdentifyDomains

IdentifyDomains = function(est, opt, full_matrix)
{
  n_matrices = dim(est$solObj)[1]
  n_rows = nrow(est$y)
  n_cols = ncol(est$y)
  dom = {}
  dom$scores = {}                              # boundary scores: all methods
  dom$bscores = matrix(0,n_rows,n_matrices)    # "normalized" boundary scores, only selected method
  dom$E = array(0,dim=c(n_rows,n_matrices))
  rownames(dom$E) = rownames(est$y)
  for (k in 1:n_matrices) {
    if (opt$verbose) write(paste('-- matrix #',k,'...',sep=''),stderr())

    # compute boundary scores
    if (full_matrix) {
      z = est$solObj[k,,]
      rownames(z) = rownames(est$y)
      colnames(z) = colnames(est$y)
    } else {
      if (opt$verbose) write('Inverse-rotating input matrices...',stderr())
      z = MatrixInverseRotate45(est$solObj[k,,])
      rownames(z) = colnames(z) = rownames(est$y)
    }
    dom$scores[[k]] = MatrixBoundaryScores(z,distance=opt$distance,d2=opt$distance2,skip=opt$'skip-distance')       # non-normalized boundary scores, all methods
    dom$scores[[k]][est$ignored_rows,] = NA                                                                         # ignored rows have no score (i.e. NA)
    if (opt$method=='inter') dom$bscores[,k] = max(dom$bscores[,k],na.rm=TRUE)-dom$bscores[,k]                      # reverse inter score
    dom$bscores[,k] = local_maxima_score(dom$scores[[k]][,opt$method],maxd=opt$'flank-dist',scale=FALSE)            # local normalization (no scaling to [0,1])
    dom$bscores[,k] = dom$bscores[,k]/quantile(dom$bscores[,k],probs=0.99,na.rm=TRUE)                               # scale to top 1% (TODO: is this necessary, or should we just replace with max?)

    # identify local maxima
    lmax = local_maxima(as.vector(dom$bscores[,k]),tolerance=opt$tolerance,alpha=opt$alpha,maxd=opt$'flank-dist')
    dom$E[,k] = lmax
  }
  rownames(dom$bscores) = rownames(est$y)
  
  return(dom)
}




# ########################
#  OPERATION = DOMAINS
# ########################

op_domains <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
    make_option(c("--min-lambda"), default=0.0, help="Minimum value for parameter lambda (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--max-lambda"), default=1.0, help="Maximum value for parameter lambda (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--n-lambda"), default=2, help="Number of lambdas (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--log2-lambda"), action="store_true",default=FALSE, help="Use log2 scale for lambda range (only applicable if estimated max-lambda was set to Inf)."),
    make_option(c("--gamma"), default=0.0, help="Value for sparsity parameter gamma (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--row-labels"), action="store_true",default=FALSE, help="Input matrix has row labels"),
    make_option(c("--ignored-loci"), default="", help="Ignored row and column names found in this list (not required) [default=\"%default\"]."),
    make_option(c("--preprocess"), default="none", help="Matrix preprocessing: none (default), max, mean, log2, log2mean, rank, dist, distlog2."),
    make_option(c("--distance"), default=5, help="Distance from diagonal (in number of bins) used for boundary score computation [default \"%default\"]."),
    make_option(c("--distance2"), default=0, help="Distance from diagonal (in number of bins) used for boundary score computation [default \"%default\"]."),
    make_option(c("--skip-distance"), default=1, help="Distance from diagonal (in number of bins) to be skipped [default \"%default\"]."),
    make_option(c("--method"), default="ratio", help="Boundary score method: ratio, diffratio, intra-max, intra-min, intra-right, intra-left, inter, diff, DI, product-max, product-min [default \"%default\"]."),
    make_option(c("--tolerance"), default=0.01, help="Percent difference cutoff for merging local maxima [default \"%default\"]."),
    make_option(c("-a","--alpha"), default=0.10, help="Minimum difference by which local maxima are greater than neighboring values [default \"%default\"]."),
    make_option(c("--flank-dist"), default=10, help="Local maxima neighborhood radius (in number of bins) [default \"%default\"]."),
    make_option(c("--track-dist"), default=10, help="Maximum distance (number of bins) from diagonal for track generation  [default=%default]."),
    make_option(c("--bins"), default="", help="Comma-separated bins to be highlighted [default \"%default\"]."),
    make_option(c("--presentation"), default="none", help="Presentation style: detailed, tracks [default \"%default\"].")
  )
  usage = 'hic-matrix.r domains [OPTIONS] MATRIX(tsv/RData)';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  if (opt$verbose) print_options(opt)
  files <- arguments$args;
  if (length(files)!=1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  f <- files[1]
  out_dir <- opt$'output-dir'
  presentation <- opt$'presentation'

  # parameters  
  bins = as.integer(strsplit(opt$'bins',',')[[1]])
  if (length(bins)==0) { bin_points = c() } else { bin_points = (bins-min(bins))/(max(bins)-min(bins)) }
  
  # create output directory
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Error: output directory already exists!',stderr()); quit(save='no'); }

  # load data
  if (opt$verbose) write("Loading data...",stderr())
  est = LoadEstimation(f,opt,replace.na=TRUE)
  if (is.null(est$lambdas)==TRUE) { write("Error: maximum lambda cannot be infinite, use the 'extract' operation to select specific lambdas.",stderr()); quit(save='no') }

  # compute boundary scores
  if (opt$verbose) { write("Computing boundary scores...",stderr()); }
  full_matrix = is_full_matrix(est$y)
  dom = IdentifyDomains(est,opt,full_matrix)
  
  # setup presentation style
  n_matrices = dim(est$solObj)[1]
  n_rows = nrow(est$x)
  n_cols = ncol(est$x)
  
  if (presentation!="none") {
    if (presentation=="tracks") {
      n_panels = 10
      a = n_rows/(n_panels*opt$'track-dist')
      pdf(paste(out_dir,'/domains.pdf',sep=''),height=9,width=9*a)
      par(mfrow=c(n_panels,1),mar=c(0.5,0.5,0.5,0.5))
    } else { 
      pdf(paste(out_dir,'/domains.pdf',sep=''))
    }
 
    # generate plots
    if (opt$verbose) { write("Generating plots...",stderr()); }
    I = as.integer(seq(1,n_rows,length.out=200))         # sample matrix for heatmap generation
    for (ll in 1:n_matrices) {
      if (opt$verbose) { write(paste('-- matrix #',ll,'...',sep=''),stderr()); }

      if (presentation=="heatmaps_only") par(mfrow=c(1,1))
      else if (presentation=="detailed") par(mfrow=c(2,2))

      # preprocess input matrices
      if (full_matrix==TRUE) { z = MatrixRotate45(est$solObj[ll,,],opt$'track-dist')
      } else { z = est$solObj[ll,,] }
      if (opt$'preprocess'=='none') z = log2(z+1)
      if (ll==1) z1 = z
    
      # plot heatmap and domains
      lmax = dom$E[,ll]
      s0 = as.vector(dom$scores[[ll]][,opt$method])     # non-normalized boundary score
      s = as.vector(dom$bscores[,ll])                   # normalized boundary score
      enorm = which(lmax==1)/length(lmax)
      if (presentation=="detailed") {
        # TODO: this needs the full matrices to work properly...
#        z1 = z1[I,I]; z1[z1<=0]=0; image(max(z1)-z1,main='Raw Hi-C heatmap')
#        rect(enorm[-length(enorm)],enorm[-length(enorm)],enorm[-1],enorm[-1],col=NULL,border='cyan',lwd=2)
#        z = z[I,I]; z[z<=0]=0; image(max(z)-z,main=substitute(paste('Estimated Hi-C matrix for ',lambda,'=',LAMBDA,sep=''),list(LAMBDA=round(est$lambdas[ll],3))))
#        rect(enorm[-length(enorm)],enorm[-length(enorm)],enorm[-1],enorm[-1],col=NULL,border='cyan',lwd=2)
#        plot(s,col='black',type='l',main='Boundary scores across chromosome',xlab='position (number of bins)',ylab='boundary score');
#        pp=s; pp[lmax==0]=NA; points(pp,col='cyan',pch=20,cex=1)
#        boxplot(s,n=10,main='Histogram of boundary scores',ylab='boundary score')
      } else if (presentation=="tracks") {
        image(z,xaxt='n',yaxt='n')
        lines(seq(0,1,length.out=length(s0)),(s0-min(s0,na.rm=TRUE))/(max(s0,na.rm=TRUE)-min(s0,na.rm=TRUE)),col='green4')
        lines(seq(0,1,length.out=length(s)),s,col='blue')
        abline(v=(which(lmax==1)-1)/(length(lmax)-1),col='blue')
        d = 1/n_rows
        for (pp in bin_points) rect(pp-d/2,0,pp+d/2,0.02,col='purple',border='purple')
      }
    
    }
    dev.off()
  }
  
  if (opt$verbose) { write("Saving data...",stderr()); }

  # create boundary files for each lambda
  for (k in 1:n_matrices) {
    f = paste(out_dir,'/boundaries.k=',formatC(k,width=3,format='d',flag='0'),'.tsv',sep='')
    w = which(dom$E[,k]==1)
    write.table(t(rbind(rownames(dom$E)[w],round(dom$bscores[w,k],6))),file=f,sep='\t',row.names=F,col.names=F,quote=F)
  }
    
  # create boundary score files (all scoring methods) for each lambda
  for (k in 1:n_matrices) {
    f = paste(out_dir,'/all_scores.k=',formatC(k,width=3,format='d',flag='0'),'.tsv',sep='')
    score_table = cbind(rownames(dom$scores[[k]]),round(dom$scores[[k]],6))
    colnames(score_table) = c("locus",colnames(dom$scores[[k]]))
    write.table(score_table,file=f,sep='\t',row.names=F,col.names=T,quote=F)
  }
    
  # store boundary score data (all lambdas in one file)    
  write.table(round(dom$bscores,6),file=paste(out_dir,'/boundary_scores.tsv',sep=''),row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

  # store extrema data (all lambdas in one file)
  write.table(dom$E,file=paste(out_dir,'/extrema.tsv',sep=''),row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

  # save data
  save(opt,est,dom,file=paste(out_dir,'/domains.RData',sep=''))
  
  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}




# ####################### #
#  scale_vector           #
# ####################### #

scale_vector = function(x,max_val)
{
  y = x
  y[y>max_val] = max_val
  y = y/max_val
  return(y)
}



# ####################### #
#  local_maxima_score     #
# ####################### #

local_maxima_score = function(x,maxd,scale)
{
  n = length(x)
  y = rep(0,n)
  for (k in 1:n) {
    if (is.na(x[k])) { y[k] = NA
    } else {
      min_left = max(k-maxd,1)+which.min(x[max(k-maxd,1):k])-1    # minimum value to the left of k
      min_right = k+which.min(x[k:min(k+maxd,length(x))])-1       # minimum value to the right of k
      y[k] = min(abs(x[k]-x[min_left]),abs(x[k]-x[min_right]))
    }
  }  
  if (scale==TRUE) y = y/max(y,na.rm=TRUE)
  return(y) 
}



# ###################
#  local_maxima     #
# ###################

local_maxima = function(x,tolerance,alpha,maxd)
{
  y = rep(0,length(x))
  x_order = order(x,decreasing=TRUE,na.last=NA)
  kk = 1
  while (kk<=length(x_order)) {
    k = x_order[kk]
    min_left = max(k-maxd,1)+which.min(x[max(k-maxd,1):k])-1    # minimum value to the left of k
    min_right = k+which.min(x[k:min(k+maxd,length(x))])-1       # minimum value to the right of k
    for (j in min_left:min_right) if (y[j]==0) y[j] = -1
    if (x[k]>=alpha) y[k] = 1
    while ((kk<=length(x_order))&&(y[x_order[kk]]!=0)) kk = kk + 1
  }  

  # add boundaries at NA values  
  d = which(is.na(x))
  di = which(d[-1]-d[-length(d)]>1)
  y[sort(c(d[di]+1,d[di+1]-1))] = 1

  #plot(x,type='l'); i=which(y==1); points(i,x[i],col='red',pch=18); j=which(y==-1); #points(j,x[j],col='blue',pch=18);  
  y[y==-1] = 0

  return(y) 
}






# ##########################
#  OPERATION = DOMAIN-DIFF
# ##########################

op_domain_diff <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
    make_option(c("--row-labels"), action="store_true",default=FALSE, help="Input matrix has row labels"),
    make_option(c("--min-lambda"), default=0.0, help="Minimum value for parameter lambda (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--max-lambda"), default=1.0, help="Maximum value for parameter lambda (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--n-lambda"), default=2, help="Number of lambdas (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--log2-lambda"), action="store_true",default=FALSE, help="Use log2 scale for lambda range (only applicable if estimated max-lambda was set to Inf)."),
    make_option(c("--gamma"), default=0.0, help="Value for sparsity parameter gamma (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--ignored-loci"), default="", help="Ignored row and column names found in this list (not required) [default=\"%default\"]."),
    make_option(c("--preprocess"), default="none", help="Matrix preprocessing: none (default), max, mean, log2, log2mean, rank, dist, distlog2."),
    make_option(c("--distance"), default=5, help="Distance from diagonal (in number of bins) used for boundary score computation [default \"%default\"]."),
    make_option(c("--distance2"), default=0, help="Distance from diagonal (in number of bins) used for boundary score computation [default \"%default\"]."),
    make_option(c("--skip-distance"), default=1, help="Distance from diagonal (in number of bins) to be skipped [default \"%default\"]."),
    make_option(c("--method"), default="ratio", help="Boundary score method: ratio, diffratio, intra-max, intra-min, intra-right, intra-left, inter, diff, DI [default \"%default\"]."),
    make_option(c("--tolerance"), default=0.01, help="Percent difference cutoff for merging local maxima [default \"%default\"]."),
    make_option(c("--alpha"), default=0.25, help="Minimum difference by which local maxima are greater than neighboring values [default \"%default\"]."),
    make_option(c("--delta"), default=0.25, help="Minimum difference for calling differential boundaries [default \"%default\"]."),
    make_option(c("--flank-dist"), default=10, help="Local maxima neighborhood radius (in number of bins) [default \"%default\"]."),
    make_option(c("--track-dist"), default=10, help="Maximum distance (number of bins) from diagonal for track generation  [default=%default]."),
    make_option(c("--bins"), default="", help="Comma-separated bins to be highlighted [default \"%default\"].")
  )
  usage = 'hic-matrix.r domain-diff [OPTIONS] MATRIX-FILES(tsv/RData)'
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  if (opt$verbose) print_options(opt)
  files <- arguments$args;
  if (length(files)<1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  out_dir <- opt$'output-dir'
  bins = as.integer(strsplit(opt$'bins',',')[[1]])
  if (length(bins)==0) { bin_points = c() } else { bin_points = (bins-min(bins))/(max(bins)-min(bins)) }
  
  # create output directory
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no') }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Error: output directory already exists!',stderr()); quit(save='no') }

	# loading matrix data
	est = {}
	dom = {}
	for (f in 1:length(files)) { 
    if (opt$verbose) write(paste('Loading input matrix ',files[f],'...',sep=''),stderr())
	  est[[f]] = LoadEstimation(files[[f]],opt,replace.na=TRUE)
	  if (is.null(est[[f]]$lambdas)==TRUE) { 
	    write("Error: maximum lambda cannot be infinite, use the 'extract' operation to select specific lambdas.",stderr())
	    quit(save='no') 
	  }
	  if (opt$verbose) { write("Computing boundary scores...",stderr()); }
	  full_matrix = is_full_matrix(est[[f]]$y)
	  dom[[f]] = IdentifyDomains(est[[f]],opt,full_matrix)
	}

  # repeat for every lambda!
  for (ll in 1:length(est[[1]]$lambdas)) {
	  if (opt$verbose) write(paste('Processing lambda=',est[[1]]$lambdas[ll],'...',sep=''),stderr())
  
		# obtaining boundary scores
		scores0 = {}
		scores = {}
		lmax = {}
		for (f in 1:length(files)) { 
		  scores0[[f]] = dom[[f]]$scores[[ll]][,opt$method]
		  scores[[f]] = dom[[f]]$bscores[,ll]
		  lmax[[f]] = dom[[f]]$E[,ll]
		}

		# process data
		if (opt$verbose) write("Identifying boundary changes...",stderr())
		d = scores[[2]]-scores[[1]]          # lmax score differences
		d2 = rep(0,length(d))
		for (i in 1:length(d2)) {
		  flank = max(1,i-opt$'flank-dist'):min(length(d2),i+opt$'flank-dist')
		  if ((lmax[[1]][i]==1)&(max(lmax[[2]][flank],na.rm=TRUE)==0)&(d[i]<=-opt$delta)) d2[i] = -1       # TODO: is this too stringent??
		  if ((lmax[[2]][i]==1)&(max(lmax[[1]][flank],na.rm=TRUE)==0)&(d[i]>=+opt$delta)) d2[i] = +1
		}
		b_gain = d2==+1
		b_loss = d2==-1
		if (opt$verbose) { write(paste('Losses = ',sum(b_loss),sep=''),stderr()); write(paste('Gains = ',sum(b_gain),sep=''),stderr()) }

		# storing results
		if (opt$verbose) { write("Storing results...",stderr()); }
		table = cbind(names(scores[[1]]),d2)
		colnames(table) = c('locus','gain/loss')
		for (f in 1:length(files)) { table = cbind(table,lmax[[f]]); colnames(table)[ncol(table)] = paste('boundary-sample-',f,sep='') }
		for (f in 1:length(files)) { table = cbind(table,round(scores[[f]],4)); colnames(table)[ncol(table)] = paste('lmax-score-sample-',f,sep='') }
		for (f in 1:length(files)) { table = cbind(table,round(scores0[[f]],4)); colnames(table)[ncol(table)] = paste('lmax-score0-sample-',f,sep='') }
		write.table(table,col.names=T,row.names=F,quote=F,sep='\t',file=paste(out_dir,'/table.k=',formatC(ll,width=3,format='d',flag='0'),'.tsv',sep=''))
		#save(opt,dom,scores,lmax,d,d2,file=paste(out_dir,'/diff.k=',formatC(ll,width=3,format='d',flag='0'),'.RData',sep=''))   # TODO: make this an option?
		
		# generate snapshots
		if (length(data)!=0) {
		  if (opt$verbose) write("Generating snapshots...",stderr())
		  pdf(paste(out_dir,'/diff.k=',formatC(ll,width=3,format='d',flag='0'),'.pdf',sep=''),height=9,width=18)
		  par(mfcol=c(6,2),mar=c(0.5,0.5,0.5,0.5))
		  flank = 2.5*opt$'track-dist'
		  s0min = min(sapply(scores0,min,na.rm=TRUE))
		  s0max = max(sapply(scores0,max,na.rm=TRUE))
		  zlim = {}
		  for (f in 1:length(files)) {
		    zlim[[f]] = max(est[[f]]$solObj[ll,,])
		    if (opt$preprocess=='none') zlim[[f]] = log2(zlim[[f]]+1)
		  }
		  for (j in which(b_loss|b_gain)) if ((j>flank)&(j<=length(d2)-flank)) {
		    if (opt$verbose) write(paste(j,lmax[[1]][j],lmax[[2]][j]),stderr())
		    I = (j-flank):(j+flank)
		    for (f in 1:length(files)) {
		      full_matrix = is_full_matrix(est[[f]]$y)
		      if (full_matrix) { mat = MatrixRotate45(est[[f]]$solObj[ll,I,I],opt$'track-dist')
		      } else { mat = est[[f]]$solObj[ll,I,1:min(ncol(est[[f]]$y),opt$'track-dist')] }
		      if (opt$preprocess=='none') mat = log2(mat+1)
		      image(mat,xaxt='n',yaxt='n',zlim=c(0,zlim[[f]]))
		      s0 = as.vector(scores0[[f]][I])
		      s = as.vector(scores[[f]][I])
		      lines(seq(0,1,length.out=length(s0)),(s0-s0min)/(s0max-s0min),col='green4')
		      lines(seq(0,1,length.out=length(s)),s,col='blue')
		      enorm = which(lmax[[f]][I]==1)/length(I)
		      abline(v=(which(lmax[[f]][I]==1)-1)/(length(I)-1),col='blue')
		      d = 1/length(I)
		      pp = 0.5
		      rect(pp-d,0,pp+d,0.04,col='purple',border='purple')
		    }
		  }
		  dev.off();
		}

  }
  		
  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no')
}





# #####################
#  OPERATION = CORR
# #####################

op_corr <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-file"), default="", help="Output text file (required) [default \"%default\"].")
  );
  usage = 'hic-matrix.r corr [OPTIONS] MATRIX-TSV-1 MATRIX-TSV-2 ...';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
  opt <- arguments$options
  files <- arguments$args
  if (length(files)<2) { write(paste('Usage:',usage),stderr()); quit(save='no') }
  fout <- opt$'output-file'

  # reading matrices
  M = {}
  for (k in 1:length(files)) { 
    if (opt$verbose) write(paste('Reading ',files[k],'...',sep=''),stderr())
    M[[files[k]]] = as.matrix(read.table(files[k]))
  }
  
  # score for each matrix row, then combine across matrices
  score = apply(sapply(M,function(m) apply(m,1,max,na.rm=TRUE)),1,min,na.rm=TRUE)

  # compute correlations for different cutoffs
  write('RESULTS',fout)
  for (cutoff in unique(quantile(score,probs=c(0.01,0.05,0.10),na.rm=TRUE))) {
    # filter
    z = score>cutoff
    write(paste('Correlations for cutoff=',cutoff,' (removed ',sum(z==FALSE),' rows)...',sep=''),fout,append=TRUE)

    # convert to vectors
    V = sapply(M,function(m) as.vector(m[z,z]))

    # correlations
    C = {}
    for (method in c('pearson','spearman')) C[[method]] = cor(V,use='complete.obs',method=method)

    # print
    for (method in names(C)) {
      write(paste('Correlation method = ',method,sep=''),file=fout,append=TRUE)
      write.table(cbind(round(C[[method]],3),colnames(M)),sep='\t',quote=F,row.names=FALSE,col.names=FALSE,file=fout,append=TRUE)
      write('',file=fout,append=TRUE)
    }
  }
  
  # done
  if (opt$verbose) write('Done.',stderr())
  quit(save='no')
}





# ########################
#  OPERATION = COMPARE
# ########################

op_compare <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-file"), default="", help="Output pdf file (required) [default \"%default\"]."),
    make_option(c("--bin-size"), default=0, help="Bin size in nucleotides for RPKB calculation [default \"%default\"]."),
    make_option(c("--bins"), default="", help="Comma-separated list of bins to be included in the snapshot [default \"%default\"]."),
    make_option(c("--min-distance"), default=0, help="Minimum distance in nucleotides [default \"%default\"]."),
    make_option(c("--min-lambda"), default=0.0, help="Minimum value for parameter lambda (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--max-lambda"), default=1.0, help="Maximum value for parameter lambda (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--n-lambda"), default=2, help="Number of lambdas (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--log2-lambda"), action="store_true",default=FALSE, help="Use log2 scale for lambda range (only applicable if estimated max-lambda was set to Inf)."),
    make_option(c("--gamma"), default=0.0, help="Value for sparsity parameter gamma (only applicable if estimated max-lambda was set to Inf) [default=%default]."),
    make_option(c("--max-dist"), default=0, help="(Only for tsv files) Maximum distance from diagonal in number of bins [default \"%default\"]."),
    make_option(c("--n-dist"), default=10, help="(Only for tsv files) Number of distances to be tested [default \"%default\"].")
  )
  usage = 'hic-matrix.r compare [OPTIONS] MATRIX-1 MATRIX-2 (* matrices can be tsv or estimated RData)'
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  fraction <- opt$"contact-map-fraction";
  files <- arguments$args;
  if (length(files)!=2) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  f1 <- files[1]
  f2 <- files[2]
  fout <- opt$'output-file'
  bin_size = opt$'bin-size'
  bins = as.integer(strsplit(opt$bins,',')[[1]])
  min_dist = opt$'min-distance'

  # Error checking on input arguments
  if (fout=="") { write('Error: please specify output file!',stderr()); quit(save='no'); }
  
  ext = sub("^.*[.]","",f1)
  if (ext=="RData") {
    # compute and plot correlations
    if (opt$verbose) { write("Loading data...",stderr()); }
    e1 = ExtractEstimations(f1,opt)
    e2 = ExtractEstimations(f2,opt)
    if (prod(e1$lambdas==e2$lambdas)!=1) { write('Error: different lambda range in the two samples!',stderr()); quit(save='no'); }
    if (prod(e1$gammas==e2$gammas)!=1) { write('Error: different gamma range in the two samples!',stderr()); quit(save='no'); }
    if (nrow(e1$y)!=nrow(e2$y)) { write('Error: different matrix sizes in samples!',stderr()); quit(save='no'); }
    if (e1$opt$preprocess!=e2$opt$preprocess) { write('Error: the two matrices have been estimated using a different preprocessing method!',stderr()); quit(save='no'); }
    if (tolower(e1$opt$algorithm)!=tolower(e2$opt$algorithm)) { write('Error: the two matrices have been estimated using a different method!',stderr()); quit(save='no'); }
    n_lambda = n_matrices = length(e1$lambdas);
    lambdas = e1$lambdas
    gammas = e1$gammas
    ignored_rows = unique(c(e1$ignored_rows,e2$ignored_rows))
    ignored_cols = unique(c(e1$ignored_cols,e2$ignored_cols))

    # determine off-diagonal distances
    max_dist = ncol(e1$y)
    if ((tolower(e1$opt$algorithm)=='fused2dzone_flsa')||(tolower(e1$opt$algorithm)=='fused1dzone_flsa')) {
      if (e1$opt$'zone-size'!=e2$opt$'zone-size') { write('Error: the two matrices have been estimated using a different zone size!',stderr()); quit(save='no'); }
      max_dist = as.integer(e1$opt$'zone-size')
    }
    a = c(1.0)   #c(0.5,1.0)
    distances = unique(as.integer(a*max_dist))

    # open output PDF file
    pdf(paste(fout,'.pdf',sep=''));
  
    # plot correlations
    if (opt$verbose) { write("Computing/plotting correlations...",stderr()); }
    par(mfrow=c(2,2));
    c_pearson = c_spearman = matrix(0,length(distances),n_matrices)
    rownames(c_pearson) = rownames(c_spearman) = distances
    c_pearson = compare_matrices(e1$solObj,e2$solObj,distances=distances,method='pearson',verbose=TRUE)
    c_spearman = compare_matrices(e1$solObj,e2$solObj,distances=distances,method='spearman',verbose=TRUE)
    xaxis_n = n_lambda;
    xaxis_values = e1$lambdas;
    xaxis_label = 'lambda';
    q = unique(as.integer(seq(1,xaxis_n,length.out=10)))
    qlab = round(xaxis_values[q],2);
    legend_pos = "bottomleft"
    plot_matrix(c_pearson,main='Pearson correlation',xlab=xaxis_label,ylab='correlation',q=q,qlab=qlab,legend_pos=legend_pos)
    plot_matrix(c_spearman,main='Spearman correlation',xlab=xaxis_label,ylab='correlation',q=q,qlab=qlab,legend_pos=legend_pos)

    # boxplots of matrix values
    L1 = as.list(data.frame(apply(e1$solObj,1,as.vector)))
    L2 = as.list(data.frame(apply(e2$solObj,1,as.vector)))
    boxplot(c(L1,L2),outline=FALSE)
  
    # plot degrees of freedom
    if (opt$verbose) { write("Computing/plotting degrees of freedom...",stderr()); }
    par(mfrow=c(2,2));
    N = sum(abs(row(e1$y)-col(e1$y))<=max_dist);
    tolerance = 0;                     # TODO: make this a parameter?
    df1 = df2 = rep(0,n_matrices)
    for (k in 1:n_matrices) {  
      if (opt$verbose) { write(paste('-- matrices #',k,sep=''),stderr()); }
      df1[k] = CalcDF(e1$solObj[k,,],distance=max_dist,tolerance=tolerance*mean(e1$solObj[k,,]))
      df2[k] = CalcDF(e2$solObj[k,,],distance=max_dist,tolerance=tolerance*mean(e2$solObj[k,,]))
    }
    plot(df1,type='l',col='red',xlab=xaxis_label,ylab='degrees of freedom',xaxt='n',ylim=c(0,max(df1,df2))); lines(df2,col='green4'); axis(1,at=q,qlab);
    legend('topright',c('matrix1','matrix2'),pch=16,col=c('red','green4'),inset=0.05);
    plot(df1/N,type='l',col='red',xlab=xaxis_label,ylab='degrees of freedom (normalized)',xaxt='n',ylim=c(0,1)); lines(df2/N,col='green4'); axis(1,at=q,qlab);
    legend('topright',c('matrix1','matrix2'),pch=16,col=c('red','green4'),inset=0.05);
  
    # find optimal lambda and print results
    if (opt$verbose) { write("Finding optimal lambda...",stderr()); }
    fdiff = function(x) { x[-1]-x[-length(x)] }
    df1_diff = fdiff(df1)
    df2_diff = fdiff(df2)
    plot(df1_diff,type='l',col='red',xlab=xaxis_label,ylab='df diff',xaxt='n',ylim=c(min(df1_diff,df2_diff),max(df1_diff,df2_diff))); lines(df2_diff,col='green4'); axis(1,at=q,qlab)
    optimal1 = 1 + order(df1_diff)[1]    # Finds maximum drop in df1
    optimal2 = 1 + order(df2_diff)[1]    # Finds maximum drop in df2
    optimal = min(optimal1,optimal2);
    write(paste('Correlations of matrices for lambda=',lambdas[1],sep=''),file=stdout())
    write(c(c_pearson[,1],c_spearman[,1]),file=stdout())
    write(paste('Correlations of matrices for optimal lambda=',lambdas[optimal],': ',sep=''),file=stdout())
    write(c(c_pearson[,optimal],c_spearman[,optimal]),file=stdout())

    # plot value of objective function
    if (length(lambdas)>1) {
      if (opt$verbose) { write("Computing penalty-vs-loss ratios...",stderr()); }
      objf1 = rep(0,length(lambdas));
      objf2 = rep(0,length(lambdas));
      for (i in 1:length(lambdas)) {  
        if (opt$verbose) { write(paste('-- matrices #',i,sep=''),stderr()); }
        objf1[i] = CalcPenaltyLossRatio(e1$y,e1$solObj[i,,],lambdas[i],gammas);
        objf2[i] = CalcPenaltyLossRatio(e2$y,e2$solObj[i,,],lambdas[i],gammas);
      }
      ymax = max(max(objf1,na.rm=TRUE),max(objf2,na.rm=TRUE));
      plot(objf1,type='l',col='red',xlab=xaxis_label,ylab='penalty-vs-loss ratio',xaxt='n',ylim=c(0,ymax)); lines(objf2,col='green4'); axis(1,at=q,qlab);
    }  
    dev.off();

    # save correlation matrices  
    if (opt$verbose) { write("Saving correlation matrices...",stderr()); }
    Cout = cbind(lambdas,t(c_pearson))
    colnames(Cout)[1] = 'lambda'
    write.table(Cout,quote=F,row.names=F,sep='\t',file=paste(fout,'.cor.pearson.tsv',sep=''));
    Cout = cbind(lambdas,t(c_spearman))
    colnames(Cout)[1] = 'lambda'
    write.table(Cout,quote=F,row.names=F,sep='\t',file=paste(fout,'.cor.spearman.tsv',sep=''));
  
    # save
    if (opt$verbose) { write("Saving data...",stderr()); }
    save(lambdas,gammas,c_pearson,c_spearman,df1,df2,objf1,objf2,optimal,file=paste(fout,'.RData',sep=''));

  } else {
    # assume matrices are in tsv format
    mat1 = as.matrix(read.table(files[1],row.names=1,check.names=F))
    mat2 = as.matrix(read.table(files[2],row.names=1,check.names=F))
    if ((nrow(mat1)!=nrow(mat2))||(ncol(mat1)!=ncol(mat2))) { write('Error: different matrix sizes in samples!',stderr()); quit(save='no'); }

    # replace NAs with 0
    mat1[is.na(mat1)] = 0
    mat2[is.na(mat2)] = 0

    # distance from diagonal
    max_dist = ncol(mat1)
    if ((opt$"max-dist">0)&&(opt$"n-dist">0)) max_dist = min(max_dist,opt$"max-dist")
    full_matrix = is_full_matrix(mat1)
    D = abs(row(mat1)-col(mat1)) + 1

    # compute correlations
    methods = c("pearson","spearman")
    distances = rev(as.integer(1:opt$"n-dist"*(max_dist/opt$"n-dist")))
    for (m in methods) {
      C = matrix(0,1,length(distances))
      colnames(C) = c(paste("d=",distances,sep=''))
      for (k in 1:length(distances)) {
        d = distances[k]
        C[k] = ifelse(full_matrix==TRUE,cor(as.vector(mat1[D<=d]),as.vector(mat2[D<=d]),method=m),cor(as.vector(mat1[,1:d]),as.vector(mat2[,1:d]),method=m))
      }
      C_table = cbind("0",round(C,4))
      colnames(C_table) = c("lambda",colnames(C))
      write.table(C_table,quote=FALSE,row.names=FALSE,sep='\t',file=paste(fout,'.cor.',m,'.tsv',sep=''))
    }
  }
  
  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}




# ########################
#  OPERATION = COMPARE2
# ########################

op_compare2 <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-file"), default="", help="Output pdf file (required) [default \"%default\"]."),
    make_option(c("-b","--bin-size"), default=0, help="Bin size in nucleotides [default \"%default\"].")
  );
  usage = 'hic-matrix.r compare2 [OPTIONS] MATRIX1 MATRIX2';
  
  # get command line options (if help option encountered print help and exit)
  arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
  opt <- arguments$options;
  files <- arguments$args;
  if (length(files)!=2) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # input arguments
  f1 <- files[1];
  f2 <- files[2];
  fout <- opt$'output-file';
  bin_size = as.integer(opt$'bin-size')
  by_dist = 4000000/bin_size
  
  # Error checking on input arguments
  if (bin_size<=0) { write('Error: please specify a valid bin size!',stderr()); quit(save='no'); }
  if (fout=="") { write('Error: please specify output file!',stderr()); quit(save='no'); }

  # Load matrices
  if (opt$verbose) { write("Loading data...",stderr()); }
  x1 = as.matrix(read.table(f1,row.names=1,check.names=F))
  x2 = as.matrix(read.table(f2,row.names=1,check.names=F))
  if (nrow(x1)!=nrow(x2)) { write('Error: different matrix sizes in samples!',stderr()); quit(save='no'); }
  n_rows = nrow(x1)
  
  # open output PDF file
  pdf(fout);
  
  # Compute and plot correlations as a function of value
  if (opt$verbose) { write("Computing/plotting correlations as a function of value...",stderr()); }
  y1 = PreprocessMatrix(x1,preprocess='dist',pseudo=1,cutoff=0)
  y2 = PreprocessMatrix(x2,preprocess='dist',pseudo=1,cutoff=0)
  x1_rpkm = x1/(sum(x1)/1000000)/(bin_size/1000)
  x2_rpkm = x2/(sum(x2)/1000000)/(bin_size/1000)
  cutoffs = 2^seq(log2(0.01),log2(10.0),length.out=10)
  r0 = sapply(cutoffs,function(c) { z1=x1_rpkm>=c; z2=x2_rpkm>=c; sum(z1*z2>0)/min(sum(z1),sum(z2))})   #sum(z1+z2>0)})
  plot(cutoffs,r0,type='l',col='magenta',lwd='2',ylim=c(0,1),xlab='RPKM cutoff',ylab='reproducibility',main='reproducibility as a function of RPKM and enrichment cutoffs')
  r = {}
  cc = c(1.0,2.0,3.0)
  clr = c('red','green4','blue')
  for (k in 1:length(cc)) {
    r[[cc[k]]] = sapply(cutoffs,function(c) { z1=(y1>=cc[k])&(x1_rpkm>=c); z2=(y2>=cc[k])&(x2_rpkm>=c); sum(z1*z2>0)/min(sum(z1),sum(z2))})    #sum(z1+z2>0)})
    lines(cutoffs,r[[cc[k]]],col=clr[k],lwd='2')
  }
  legend("bottomright",legend=c(0,1,2,3),lty=rep(1,4),col=c('magenta',clr),text.col="black",bg="gray90")

  # Compute and plot correlations as a function of distance
  if (opt$verbose) { write("Computing/plotting correlations as a function of distance...",stderr()); }
  D = row(x1)-col(x1)
  distances = seq(0,n_rows/2,by=by_dist)
  prep = c('none','log2','distmax')
  c = {}
  for (method in c('pearson','spearman')) {
    c[[method]] = matrix(0,length(prep),length(distances)-1)
    colnames(c[[method]]) = distances[-1]
    rownames(c[[method]]) = prep
    for (j in 1:length(prep)) {
      y1 = PreprocessMatrix(x1_rpkm,preprocess=prep[j],pseudo=0.001,cutoff=0.05)
      y2 = PreprocessMatrix(x2_rpkm,preprocess=prep[j],pseudo=0.001,cutoff=0.05)
      for (k in 2:length(distances)) {
        J = (D>distances[k-1])&(D<=distances[k])
        c[[method]][j,k-1] = cor(as.vector(y1[J]),as.vector(y2[J]),method=method)
      }
    }
    plot_matrix(c[[method]],main=paste(method,' correlation',sep=''),xlab='distance',ylab='correlation',q=1:length(distances[-1]),qlab=distances[-1],legend_pos="topright")
  }
  
  # save
  if (opt$verbose) { write("Saving data...",stderr()); }
  save(cutoffs,r0,r,c,file=paste(fout,'.RData',sep=''));

  # done
  if (opt$verbose) { write("Done.",stderr()); }
  quit(save='no');
}




# ##############################
#  OPERATION = DOMAIN-CMP
# ##############################

op_domain_cmp <- function(cmdline_args) 
{
  # process command-line arguments
  option_list <- list(
    make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
    make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"].")
  )
  usage = 'hic-matrix.r domain-diff [OPTIONS] BOUNDARY-SCORE-FILES'
  
  # get command line options (if help option encountered print help and exit)
  arguments = parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
  opt = arguments$options
  files = arguments$args
  if (length(files)<1) { write(paste('Usage:',usage),stderr()); quit(save='no'); }

  # create output directory
  out_dir = opt$'output-dir'
  if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
  if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Error: output directory already exists!',stderr()); quit(save='no'); }
  
  # load data
  if (opt$verbose) { write("Loading data...",stderr()); }
  scores = {}
  for (f in 1:length(files)) scores[[f]] = as.matrix(read.table(files[f],check.names=F,header=T))

  # process data
  if (opt$verbose) write("Processing data...",stderr())
  func_eval = function(method) {
    if (opt$verbose) write(paste('-- method = ',method,sep=''),stderr())
    eval = {}
    k = 1
    for (maxd in seq(10,20,by=5)) {
      alpha = 0.1
      bscore = {}
      lmax = {}
      lmax_order = {}
      for (f in 1:length(files)) {
        bscore[[f]] = local_maxima_score(as.vector(scores[[f]][,method]),maxd=maxd,scale=FALSE)
        bscore[[f]] = bscore[[f]]/quantile(bscore[[f]],probs=0.99,na.rm=TRUE)
        lmax[[f]] = local_maxima(bscore[[f]],tolerance=0,alpha=alpha,maxd=maxd)
        lmax_order[[f]] = order(lmax[[f]]*bscore[[f]],decreasing=TRUE)
      }
      for (nbest in seq(1000,5000,by=1000)) {
        if ((sum(lmax[[1]])>=nbest)&&(sum(lmax[[2]])>=nbest)) {
          l1 = l2 = rep(0,length(lmax[[1]]))
          l1[lmax_order[[1]][1:nbest]] = 1
          l2[lmax_order[[2]][1:nbest]] = 1
		      m = sum(l1*l2)
		      u = which(l1==1)
		      u = unique(c(u-1,u,u+1))
		      u = u[(u>=1)&(u<=length(l1))]
		      m2 = sum(l2[u])
		      eval[[k]] = c( maxd, nbest, m, m2, m/nbest, m2/nbest )
		      k = k + 1
		    }
      }
    }
    eval = round(t(matrix(unlist(eval),nrow=length(eval[[1]]))),3)
    rownames(eval) = rep(method,nrow(eval))
    return(eval)
  }
  e = do.call('rbind',lapply(colnames(scores[[1]]),func_eval))
  write.table(e,quote=F,row.names=T,col.names=F,sep='\t',file=paste(out_dir,'/diff.tsv',sep=''))
  save(files,scores,e,local_maxima,file=paste(out_dir,'/diff.RData',sep=''))
  if (opt$verbose) write("Done.",stderr())
  quit(save='no')
}


  


# ########################
#  MAIN
# ########################

# process command-line arguments
args <- commandArgs(trailingOnly=T)
if (length(args)<1) {
  cat('\n');
  cat('USAGE:\n');
  cat('  hic-matrix.r OPERATION [OPTIONS] FILE-LIST\n');
  cat('\n');
  cat('VERSION:\n');
  cat(paste('  ',VERSION,'\n',sep=''));
  cat('\n');
  cat('DESCRIPTION:\n');
  cat('  Performs analysis of Hi-C sequencing contact matrix data.\n');
  cat('\n');
  cat('OPERATIONS:\n');
  cat('  preprocess   Preprocess an input matrix.\n');
  cat('  normalize    Normalize matrix by fragment effective length.\n');
  cat('  standardize  Standardize input matrix given average count statistics as a function of distance.\n');
  cat('  submatrix    Extracts submatrix from an input matrix.\n');
  cat('  stats        Computes statistics on input matrix (tsv) or estimated matrices (RData).\n');
  cat('  estimate     Estimates true contact matrix.\n');
  cat('  extract      Extracts estimated matrices for specific lambdas and gammas.\n');
  cat('  matrices     Creates matrices in text format using the estimated RData file.\n');
  cat('  corr         Computes all pairwise correlations of input matrices (tsv format).\n');
  cat('  compare      Compares estimated contact matrices.\n');
  cat('  compare2     Compares estimated contact matrices (alternative version).\n');
  cat('  heatmaps     Creates heatmaps for estimated matrices in RData file.\n');
  cat('  snapshots    Creates snapshots (triangle format) for estimated matrices in RData file(s).\n');
  cat('  transloc     Computes translocation scores.\n');
  cat('  domains      Identifies domains on input matrix (tsv) or estimated matrices (RData).\n');
  cat('  domain-diff  Identifies domain boundary differences in input samples.\n');
  cat('  domain-cmp   Compares domain boundary scores.\n');
  cat('  loops        Identify loops (i.e. interactions).\n');
  cat('  loop-diff    Computes fold-changes between corresponding elements of two matrices.\n');
  cat('\n');
  quit(save="no");
}

write(paste('Version ',VERSION,'\n',sep=''),stderr());

# install packages
write('Loading libraries...',stderr());
for (p in c("flsa","genlasso","ggplot2","optparse","pastecs","plotrix","reshape2","zoo")) 
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }

op = args[1];
args = args[seq(2,length.out=length(args)-1)];

# load libraries
ARGS <- commandArgs(trailingOnly=F);  
scriptPath <- normalizePath(dirname(sub("^--file=", "", ARGS[grep("^--file=", ARGS)])));
dyn.load(paste(scriptPath,"/hic-matrix.so",sep=''));

if (op=="preprocess") {
  op_preprocess(args);
} else if (op=="normalize") {
  op_normalize(args);
} else if (op=="standardize") {
  op_standardize(args);
} else if (op=="submatrix") {
  op_submatrix(args);
} else if (op=="transloc") {
  op_transloc(args);
} else if (op=="stats") {
  op_stats(args);
} else if (op=="estimate") { 
  op_estimate(args)
} else if (op=="extract") { 
  op_extract(args)
} else if (op=="matrices") {
  op_matrices(args);
} else if (op=="heatmaps") {
  op_heatmaps(args);
} else if (op=="snapshots") {
  op_snapshots(args);
} else if (op=="corr") {
  op_corr(args);
} else if (op=="compare") {
  op_compare(args);
} else if (op=="compare2") {
  op_compare2(args);
} else if (op=="domains") {
  op_domains(args);
} else if (op=="domain-diff") {
  op_domain_diff(args);
} else if (op=="domain-cmp") {
  op_domain_cmp(args);
} else if (op=="loops") {
  op_loops(args);
} else if (op=="loop-diff") {
  op_loopdiff(args);
} else {
  write(paste('Error: unknown operation "',op,'"!',sep=''),stderr());
}






