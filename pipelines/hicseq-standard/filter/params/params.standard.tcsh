#!/bin/tcsh

source ./inputs/params/params.tcsh

set release = inputs/release
set filter_params = "--mapq 30 --min-dist 25000 --max-offset 500 --filter-dups"

