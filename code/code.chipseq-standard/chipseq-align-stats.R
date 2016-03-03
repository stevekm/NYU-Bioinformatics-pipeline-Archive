#!/usr/bin/Rscript

## USAGE: Rscript ./code/chipseq-align-stats.R $outdir $branch "$objects"
## DESCRIPTION: create summary tables and dual barplots to visualize alignment reads
## 

# get the script arguments
args <- commandArgs(TRUE)


cat("R: dput args:") # for development & troubleshooting
cat("\n")
cat(dput(args))
cat("\n\n")

OutDir<-args[1]
Branch<-args[2]
Objects<-unlist(strsplit(args[3]," "))

cat("Outdir is ",OutDir,"",sep = "\n")
cat("Branch is ",Branch,"",sep = "\n")
cat("Objects is ","",sep = "\n")
Objects

# preallocate vectors to hold each column we will need in the resulting stats table
Total_reads <- numeric(length(Objects))
Aligned_reads <- numeric(length(Objects))
De_dup_aligns <- numeric(length(Objects))
Sample_Name <- character(length(Objects))


#
##
###
# read in the stats for each alignment
for(i in 1:length(Objects)){
  # get the sample's stats sheet
  SampleStatsSheet<-paste(Branch,as.character(Objects[i]),"stats.tsv",sep = "/")

  # add name of the sample to vector
  Sample_Name[i]<-as.character(Objects[i])
  Total_reads[i]<-read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1)["Total reads",][1]
  Aligned_reads[i]<-read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1)["Aligned reads",][1]
  De_dup_aligns[i]<-read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1)["De-duplicated alignments",][1]
  tmp_rownames<-row.names(read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1))
}


AlignmentStats<-data.frame(Total_reads,Aligned_reads,De_dup_aligns,row.names = Sample_Name)
colnames(AlignmentStats)<-gsub(pattern = " ",replacement = ".",x = tmp_rownames) # no whitespace allowed in colnames!!

cat("Alignment stats data frame is:",sep = "\n\n")
AlignmentStats

#
##
###
# Calculate the percent alignment, etc.
AlignmentStats[['Percent.Aligned.Reads']]<-signif(c(AlignmentStats[['Aligned.reads']] / AlignmentStats[['Total.reads']])*100,digits = 4 )
# calculate the percent deduplicated
AlignmentStats[['Percent.De-dup.Reads']]<-signif(c(AlignmentStats[['De-duplicated.alignments']] / AlignmentStats[['Total.reads']])*100,digits = 4 )
# calculate the number of unaligned reads
AlignmentStats[['Unaligned.Reads']]<-c(AlignmentStats[['Total.reads']] - AlignmentStats[['Aligned.reads']])
#calculate the percent unaligned reads
AlignmentStats[['Pcnt.Unaligned.Reads']]<-signif(c(AlignmentStats[['Unaligned.Reads']] / AlignmentStats[['Total.reads']])*100,digits=4)
# calculate the number of duplicated reads
AlignmentStats[['Duplicated']]<-c(AlignmentStats[['Aligned.reads']] - AlignmentStats[['De-duplicated.alignments']])
# calculate the percent of duplicated
AlignmentStats[['Percent.Dup']]<-signif(c(AlignmentStats[['Duplicated']] / AlignmentStats[['Total.reads']])*100,digits = 4 )



#
##
###
# save the values to be plotted into a transposed matrix, since thats what the barplot() likes
# # first just get the columns we want
Dup_Raw_Reads_df<-AlignmentStats[,which(colnames(AlignmentStats) %in% c("De-duplicated.alignments","Duplicated","Unaligned.Reads")) ] 
# reorder the columns because R is dumb
Dup_Raw_Reads_df<-Dup_Raw_Reads_df[c("De-duplicated.alignments","Duplicated","Unaligned.Reads")]

Dup_Raw_Reads_Matrix<-t(as.matrix(Dup_Raw_Reads_df))
# # divid the number of reads by 1million
Dup_Raw_Reads_Matrix<-signif(Dup_Raw_Reads_Matrix/1000000,digits = 4)

# # first just get the columns we want
Dup_Pcnt_Reads_df<-AlignmentStats[,which(colnames(AlignmentStats) %in% c("Percent.De-dup.Reads","Percent.Dup","Pcnt.Unaligned.Reads")) ]
# reorder
Dup_Pcnt_Reads_df<-Dup_Pcnt_Reads_df[c("Percent.De-dup.Reads","Percent.Dup","Pcnt.Unaligned.Reads")]
Dup_Pcnt_Reads_Matrix<-t(as.matrix(Dup_Pcnt_Reads_df))

#
##
###
# Set up the plots
BARPLOT_COLORS<-c("blue","purple","red")

# setup the matrix for the plot layout
Raw_Reads_Matrix_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
                                     2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                     3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
                                   .Dim = c(11L,4L), 
                                   .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))


# calculate some values for the figure margins
# based on the nchar() of the longest rowname, divided by a value
mar_divisor<-2.5
mar_widthLeft<-signif(
  max(4.1,
      max(nchar(row.names(AlignmentStats)))/mar_divisor),
  4) # 65 is too much # this may not work in RStudio but should work in pdf() or command line R 


# calculate value for the plot label scaling factor
# Names_scale<-min(0.7, 
#                max(nchar(row.names(AlignmentStats)))*.0075) # works alright up to 88 char's samplenames # this doesn't work as well, needs more tweaking
# ^ this also causes tiny labels for short names, too small, need both max and min cutoffs??
Names_scale<-0.7
# scaling factor for space between bars
Space_scale<-max(0.2, # default setting
                 nrow(AlignmentStats)*.01) # needs to work with up to 61


cat("The longest sample name is ",max(nchar(row.names(AlignmentStats))),sep = "\n")
cat("The number of samples is ",nrow(AlignmentStats),sep = "\n")


cat("mar_widthLeft is ",mar_widthLeft,"",sep = "\n")
cat("Names_scale is ",Names_scale,"",sep = "\n")
cat("Space_scale is ",Space_scale,sep = "\n")

#
##
###
# write a PDF of the plot
# pdf(file = paste0(OutDir,"/alignment_barplots",mar_divisor,"-",mar_widthLeft,".pdf"),width = 8,height = 8) # ORIGINAL
pdf(file = paste0(OutDir,"/alignment_barplots.pdf"),width = 8,height = 9)

# setup the panel layout
layout(Raw_Reads_Matrix_matrix) 
# need to set this for some reason
par(mar=c(0,0,4,0))
# call blank plot to fill the first panel
plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2) 
# set up the Legend in the first panel
legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0)
# plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
# par(mar=c(6,max(4.1,max(nchar(row.names(AlignmentStats)))/1.5),0,3)+ 0.1) # ORIGINAL 
par(mar=c(6,mar_widthLeft,0,3)+ 0.1) 

# create barplot for the two matrices
# barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Number of reads (millions)") 
# barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Percent of reads")
barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=Names_scale,xlab="Number of reads (millions)",space=Space_scale) 
barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=Names_scale,xlab="Percent of reads",space=Space_scale)

dev.off()

#
##
###
# write a CSV of the final table
# # peel off the rownames into a separate vector
SampleName<-row.names(AlignmentStats)
write.csv(x = cbind(SampleName,AlignmentStats), file = paste0(OutDir,"/alignment_stats_extended.csv"),quote = F,row.names = F)
