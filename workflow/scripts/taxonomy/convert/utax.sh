#!/usr/bin/env bash

set -euo pipefail


if [ -$# -lt 2 ]; then
    echo "usage: $0 <infile> <outfile>" >&2
    echo "  infile: *.fasta.zst" >&2
    echo "  outfile: *.fasta.zst" >&2
    exit 1
fi

infile="$1" && shift
outfile="$1" && shift

# using seqtool
zstd -dc < "$infile" |  # note '<' (necessary with symlinks)
  st replace -dr '[,:]' '_' |  # make sure that there are no ':' since they are used in the UTAX format
  st replace -dr '^\s*([a-z])__' 'tax=$1:' |  # prepare lineage: set top rank
  st replace -dr '\s*;\s*(\w)__' ',$1:' |  # replace lower ranks (no spaces should be left)
  st replace -dr '[a-z]:(,|$)' '' | st replace -dr ',$' '' |  # delete empty/unknown ranks
  st set -i '{id};{desc}' |  # append to id
  st del -d |  # delete description
  zstd -c > "$outfile"

#unique_tax.py -c tax_conflicts.fa -
