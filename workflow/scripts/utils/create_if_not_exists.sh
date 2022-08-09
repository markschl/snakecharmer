#!/usr/bin/env bash

set -e

# creates empty ZSTD archives if they don't exist

for f in "$@"; do
  if [ ! -f "$f" ]; then
    echo -n | zstd -qc > "$f"
    echo "Empty file created: $f" 1>&2
  fi
done
