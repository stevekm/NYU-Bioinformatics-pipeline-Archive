#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc/4.7.0
module load python/2.7.3

set cutoff = 0.05           # rows are ranked by average, and a fraction of the lowest scoring ones are set to zero


