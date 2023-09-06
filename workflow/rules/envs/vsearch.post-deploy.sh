#!/bin/bash

set -e

# install seqtool
"$PIPELINE_DIR/workflow/scripts/setup/seqtool.sh" $CONDA_PREFIX/bin
