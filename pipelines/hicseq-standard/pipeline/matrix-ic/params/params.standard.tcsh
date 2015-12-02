#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc
module unload python
module load python/2.7.3

set cutoff = 0.05

