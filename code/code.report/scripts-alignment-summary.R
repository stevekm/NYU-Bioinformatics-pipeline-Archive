#!/usr/bin/Rscript

## USAGE: alignment_summary_stats.R /path/to/outputdir /path/to/alignment/.../db.Rdata /path/to/projdir/aligndir
## DESCRIPTION: create summary tables and dual barplots to visualize alignment reads

# get the script arguments
args <- commandArgs()

OutDir<-args[6]
alignDB_filepath<-args[7]
AlignDir<-args[8]

# load the db.Rdata file
# # create a new environment for it
tmp_envir<-new.env()
# load the data into it
load(file =alignDB_filepath,envir =tmp_envir)
# ls(tmp_envir) # check its contents
obj_db<-tmp_envir$obj_db # # obj_db is the table with informaiton we need

# preallocate vectors to hold each column we will need in the resulting stats table; this is fastest method
Total_reads <- numeric(nrow(obj_db))
Aligned_reads <- numeric(nrow(obj_db))
De_dup_aligns <- numeric(nrow(obj_db))
Sample_Name <- character(nrow(obj_db))

# read in the stats for each alignment
for(i in 1:nrow(obj_db)){
  # sample's outdir location
  SampleAlignOutdir<-as.character(obj_db[["out-dir"]][i])
  # sample's stats sheet
  SampleStatsSheet<-paste(AlignDir,SampleAlignOutdir,"stats.tsv",sep = "/")
  # add name of the sample to vector
  Sample_Name[i]<-as.character(obj_db[["out-object"]][i])
  
  # add the values of the sample; ................................these names need to exist! ~~VVVVVVV
  Total_reads[i]<-read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1)["Total reads",][1]
  Aligned_reads[i]<-read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1)["Aligned reads",][1]
  De_dup_aligns[i]<-read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1)["De-duplicated alignments",][1]
  tmp_rownames<-row.names(read.table(file = SampleStatsSheet,header = F,sep = "\t",row.names = 1))
}

AlignmentStats<-data.frame(Total_reads,Aligned_reads,De_dup_aligns,row.names = Sample_Name)
colnames(AlignmentStats)<-gsub(pattern = " ",replacement = ".",x = tmp_rownames) # no whitespace allowed in colnames!!


# calculate the percent alignment
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

# #####################
# # Set up the plots
# 
BARPLOT_COLORS<-c("blue","purple","red")

# setup the matrix for the plot layout
Raw_Reads_Matrix_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 
                                     2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 
                                     3L, 3L, 3L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L),
                                   .Dim = c(11L,4L), 
                                   .Dimnames = list(NULL, c("V1", "V2", "V3", "V4")))

# write a PDF of the plot
pdf(file = paste0(OutDir,"/alignment_barplots.pdf"),width = 8,height = 8)
layout(Raw_Reads_Matrix_matrix) # setup the panel layout
par(mar=c(0,0,4,0)) # need to set this for some reason
plot(1,type='n',axes=FALSE,xlab="",ylab="",main = "Sequencing Reads",cex.main=2) # call blank plot to fill the first panel
legend("bottom",legend=c("Deduplicated","Duplicated","Unaligned"),fill=BARPLOT_COLORS,bty = "n",ncol=length(BARPLOT_COLORS),cex=1.0) # set up the Legend in the first panel
par(mar=c(6,max(4.1,max(nchar(row.names(AlignmentStats)))/2.5),0,3)) # plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
barplot(Dup_Raw_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Number of reads (millions)") # create barplot for the two matrices
barplot(Dup_Pcnt_Reads_Matrix,horiz = T,col=BARPLOT_COLORS,border=NA,las=1,cex.names=0.7,xlab="Percent of reads")
dev.off()

# write a CSV of the final table
# # peel off the rownames into a separate vector
SampleName<-row.names(AlignmentStats)
write.csv(x = cbind(SampleName,AlignmentStats), file = paste0(OutDir,"/alignment_stats_extended.csv"),quote = F,row.names = F)
