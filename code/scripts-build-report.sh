#!/bin/bash
# this script will generate an RMarkdown report based on the analysis from the pipeline
# currently pulling information from the ChIP-Seq pipeline
# this can be adjusted in future version
# run this script from the analysis directory; 
# # ~/projects/Collaborator_lab-ChIPseq-2015-10-09$ ./code/build_report.sh


# # # # 
# SETUP
# # # # 


# the location of the final report to output
REPORT_OUTDIR="./chipseq-analysis/report/report_output/"
FINAL_REPORT="./chipseq-analysis/report/report_output/final_report.Rmd"
touch $REPORT_OUTDIR
touch $FINAL_REPORT

# check to make sure the output report exists
if [ -f $FINAL_REPORT ];
then
	# if it already exists, delete it, and copy over the new report first page
	rm $FINAL_REPORT && cp ./report/report_building_blocks/first_page.Rmd $FINAL_REPORT
else
	# if the report does not already exist, copy over the first page of the report
	cp ./report/report_building_blocks/first_page.Rmd $FINAL_REPORT
fi


# # # # 
# BUILD REPORT
# # # # 


# for debug purposes, force the inclusion of these pages:
cat ./report/report_building_blocks/chunk_peak_calling.Rmd >> $FINAL_REPORT

# # look for alignment results : DEPRICATED
if [ -d ./demo_dirs/chipseq_standard/alignments ]
then
        # if the alignment results are present, then concatenate them to the final report
        cat ./report_building_blocks/chunk_alignment.Rmd >> $FINAL_REPORT
fi

# # look for analysis results, such as the pca results : DEPRICATED
if [ -d ./demo_dirs/chipseq_standard/pca ]
then
	# if the pca results are present, then concatenate the pca report chunk to the final report
	cat ./report_building_blocks/chunk_pca.Rmd >> $FINAL_REPORT
fi

# # look for other results, such as the call peaks data (doesn't exist!) : DEPRICATED
if [ -d ./demo_dirs/chipseq_standard/call_peaks ]
then
	cat ./report_building_blocks/chunk_peak_calling.Rmd >> $FINAL_REPORT
fi


# # # # 
# FINISH
# # # # 

# add the last page of the report
cat ./report/report_building_blocks/last_page.Rmd >> $FINAL_REPORT

# compile the report
./code/R_render_report.R $FINAL_REPORT
