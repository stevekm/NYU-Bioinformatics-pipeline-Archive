#Use Poisson regression remove systematic biases in Hi-C cis contact maps
#Ming Hu (minghu@fas.harvard.edu)
#Modified by Harris A. Lazaris (lazaris@nyu.edu)
#Last update: 08.05.2012
#Update: 12.21.2015

# Get required libraries
library(plyr)

# This script gets two arguments
# the matrix and the corresponding 
# genomic features
args <- commandArgs(trailingOnly=TRUE)

# Check number of arguments and exit
# if not correct
if (length(args) != 3) {print("USAGE: Rscript hicnorm_cis.r INPUT-MATRIX GENOMIC-FEATURES OUTPUT-MATRIX"); quit(save="no")}

#Arguments
input_matrix     <- args[1]
genomic_features <- args[2]
norm_matrix      <- args[3]

#read in input file
u<-read.table(sprintf("%s", input_matrix), header=TRUE, row.names=1, check.names=FALSE)      #user can change the name o
v<-read.table(sprintf("%s", genomic_features), header=FALSE,row.names=1)                     #user can change the name o
colnames(v) <- c("len","gcc","map")

# Create ids to merge
u1 <- u
u1$id <- rownames(u)
v1 <- v
v1$id <- rownames(v)
total <- join(u1,v1,by="id",type="left")

#change matrix into vector
u<-as.matrix(u)
# Change the diagonal to zero
diag(u) <- 0
u_vec<-u[upper.tri(u,diag=F)]

#get cov matrix
len_m<-as.matrix(log(total$len%o%total$len))
gcc_m<-as.matrix(log(total$gcc%o%total$gcc))
map_m<-as.matrix(log(total$map%o%total$map))

#Remove infinite values
len_m[!is.finite(len_m)] <- NA
gcc_m[!is.finite(gcc_m)] <- NA
map_m[!is.finite(map_m)] <- NA

#centralize cov matrix of enz, gcc
len_m<-(len_m-mean(c(len_m), na.rm=TRUE))/sd(c(len_m), na.rm=TRUE)
gcc_m<-(gcc_m-mean(c(gcc_m), na.rm=TRUE))/sd(c(gcc_m), na.rm=TRUE)

#change matrix into vector
len_vec<-len_m[upper.tri(len_m,diag=F)]
gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
map_vec<-map_m[upper.tri(map_m,diag=F)]

#fit Poisson regression: u~len+gcc+offset(map)
fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#user can use the following two lines to fit negative binomial regression: u~len+gcc+offset(map).
#library("MASS")
#fit<-glm.nb(u_vec~len_vec+gcc_vec+offset(map_vec))

#summary(fit)
coeff<-round(fit$coeff,4)
res<- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)

# TODO: make this an option 
res[is.na(res)] <- 0

#output normalized cis contact map, user can change the name of this output file
write.table(res, file=sprintf("%s", norm_matrix), row.names=TRUE, col.names=TRUE, sep="\t", quote=F)
