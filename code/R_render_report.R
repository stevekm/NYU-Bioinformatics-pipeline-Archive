#!/usr/bin/Rscript

# intercept the arguments to the script
args <- commandArgs(trailingOnly=TRUE) 

# R -e "rmarkdown::render('knitr_example.Rmd’)"
print(args)
rmarkdown::render(args)