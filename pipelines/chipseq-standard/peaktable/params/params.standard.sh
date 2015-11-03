#!/bin/bash

peakdist=500      # max distances between peaks that will be merged

annot_params="annotate2 -i --upstream-dist 100000 --downstream-dist 100000 --proximal-dist 1000 inputs/release/gene-name.bed"

