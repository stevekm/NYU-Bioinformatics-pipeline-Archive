#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload r
module load r/3.2.3

set diffbind_factor = 'group'
set diffbind_blocking_factor = ''


