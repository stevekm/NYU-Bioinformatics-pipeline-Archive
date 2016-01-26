#!/bin/Rscript

#Rscript create-topdom-matrix.r INPUT-MATRIX OUTPUT-MATRIX

#Load the required packages - install if not present
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr", repos='http://cran.us.r-project.org')}
library("stringr")
options(scipen=999)

###### FUNCTIONS ######

split <- function(x){unlist(str_split(x,':|-'))}

####################### 

#Get the arguments
args <- commandArgs(trailingOnly=TRUE)

#Check the arguments
if(length(args)!=2){print("Please provide right number of arguments."); quit(save="no")}

input_matrix <- args[1]
out_matrix   <- args[2]

df1 <- read.table(sprintf("%s", input_matrix), header=TRUE, row.names=1, stringsAsFactors=FALSE)
rows <- rownames(df1)

#Prepare vectors to save the components of rhe rownames
chroms <- c()
start  <- c()
end    <- c()

for (i in 1:length(rows)){chroms[i] <- print(split(rows[i])[1]); start[i] <- print(split(rows[i])[2]); end[i] <- print(split(rows[i])[3])}
df2 <- data.frame(cbind(chroms,as.numeric(start),as.numeric(end),df1))
df2[,2] <- df2[,2] - 1

# Write the output
write.table(df2,out_matrix,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
