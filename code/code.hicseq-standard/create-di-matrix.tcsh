#!/bin/tcsh

##
## USAGE: create-di-matrix.tcsh INPUT-MATRIX OUTPUT-MATRIX
##

set input = $1
set output = $2

cat $input | sed '1d' | sed 's/:/	/' | sed 's/-/	/' | tr ' ' '	' >! $output 
