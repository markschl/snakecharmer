#!/usr/bin/env bash

# Downloads a reference database from a specified location and creates a compressed
# FASTA file containing sequence IDs and QIIME-formatted taxonomy as description line.
# The taxonomy will have the whole lineage defined (no ranks missing). If names
# are not known, then they should be empty.

set -euo pipefail


if [ $# -lt 3 ]; then
    echo "usage: $0 <type> <url> <outfile>" >&2
    echo "  type: unite_otus, ..." >&2
    echo "  url: file URL" >&2
    echo "  outfile: *.fasta.zst" >&2
    exit 1
fi

type="$1" && shift
url="$1" && shift
outfile="$1" && shift


if [[ "$type" =~ ^(unite_otus)$ ]]; then
  # run sub-script
  script_dir=$(dirname $0)
  $script_dir/obtain/$type.sh "$url" "$outfile"

else
  echo "Unknown taxonomy database type: $type" >&2
  exit 1
fi
