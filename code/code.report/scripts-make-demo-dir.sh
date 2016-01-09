#!/bin/bash

rsync -a --include '*/' --exclude '*' ~/pipelines/aris_pipelines/chipseq_standard ~/pipelines/reporting/demo_dirs
