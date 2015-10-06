#!/bin/bash

# this script is not for invocation, only for reference

# use rsync to mirror directories and their contents

rsync -avK --dry-run pipelines_0/reporting/ pipelines/report/
# -a = archive mode, preserver metadata
# -v = verbose
# -K = preserves symlinks
# --dry-run = test run output, doesn't make any changes; remove this to run it for real
