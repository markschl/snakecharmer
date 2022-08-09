#!/bin/bash

set -e

# usearch
if ! command -v usearch > /dev/null; then
    echo 'USEARCH not found. Download from https://www.drive5.com/usearch, rename it to simple "usearch" and make sure it is in $PATH'
    exit 1
fi

env >&2

# seqtool
$PIPELINE_DIR/workflow/scripts/setup/seqtool.sh $CONDA_PREFIX/bin
