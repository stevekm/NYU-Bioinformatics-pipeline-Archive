#!/usr/bin/Rscript
# this script is called by build_report.sh
# this script takes an input argument in the form of a .Rmd file and 
# compiles it based on the file's YAML header

# intercept the arguments to the script
args <- commandArgs(trailingOnly=TRUE) 

print(args) # probably do not need this; outputs the input file name

# this command creates the output file
rmarkdown::render(args)


# this is an old version of the command, use this as a standalone on terminal outside of R
# R -e "rmarkdown::render('knitr_example.Rmd’)"
