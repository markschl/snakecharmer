#!/usr/bin/env bash

set -xeuo pipefail


if [ $# -lt 2 ]; then
    echo "usage: $0 <url> <outfile>" >&2
    echo "  url: file URL" >&2
    echo "  outfile: *.fasta.zst" >&2
    exit 1
fi

url="$1" && shift
outfile="$1" && shift

tmp_dir="${outfile%.*.*}"_tmp
mkdir -p "$tmp_dir"

# obtain file
echo "Downloading $url..." >&2
wget -q --no-check-certificate -O "$tmp_dir/unite.tar.gz" "$url"
mkdir -p "$tmp_dir/unite"
echo "Uncompressing..." >&2
tar -C "$tmp_dir/unite" -xzf "$tmp_dir/unite.tar.gz"

# select the dynamic clusters
# TODO: make this configurable
# TODO: allow spaces in path?
fa=$(echo "$tmp_dir/unite"/*/*_dynamic_*.fasta)
tax=$(echo "$tmp_dir/unite"/*/*_dynamic_*.txt)

echo "Selected $fa" >&2

# Add taxonomy to FASTA and do some basic cleaning
# (but no filtering) 
st set -l "$tax" -d '{l:2}' "$fa" |
  iconv -c -f UTF-8 -t ASCII |  # replace invalid ASCII
  st replace -dr ';s__[^_]+_sp\.?\s*$' ';s__unidentified' | # Convert "Genus_sp" or "Genus_sp." to unidentified (there are a few of these)
  st replace -dr '([a-z]__)unidentified' '$1' | # convert "unidentified" to empty string
  st replace -d ';' '; ' |  # TODO: remove
  zstd -c > "$outfile"

# clean up
rm -rf "$tmp_dir"
