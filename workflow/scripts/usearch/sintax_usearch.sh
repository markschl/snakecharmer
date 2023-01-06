#!/usr/bin/env bash

set -xeuo pipefail


if [ $# -lt 3 ]; then
  echo "Usage: $0 <program> <infile> <refdb> <outfile> [sintax_args]" >&2
  exit 1
fi

program="$1" && shift
infile="$1" && shift
refdb_compr="$1" && shift
outfile="$1" && shift

if ! [[ "$program" =~ ^(usearch|vsearch)$ ]]; then
  echo "SINTAX program argument needs to be either 'usearch' or 'vsearch'" >&2
  exit 1
fi

out_name=$(basename $outfile)
outdir=$(dirname $outfile)
sintax_out=$outdir/sintax/"${out_name%.*}"
refdb="${outfile%.gz}"_ref.fasta

mkdir -p $outdir/sintax

# first, uncompress the database
zstd -dqf "$refdb_compr" -o "$refdb"

$program -sintax "$infile" -db "$refdb" \
  -tabbedout "$sintax_out" \
  -quiet \
  -log /dev/stdout \
  "$@" 1>&2

# remove the uncompressed database
rm "$refdb"

# convert to Qiime format
# TODO: include correct confidence (confidence = 1 currently)
printf "Feature ID\tTaxon\tConfidence\n" | gzip -c > "$outfile"
awk -v IFS="\t" -v OFS="\t" '{print $1, $4, 1}' "$sintax_out" |
  sed -E 's/([a-z]):([^,\t]+),/\1__\2;/g' |  # convert lineages to QIIME format
  sed -E 's/([a-z]):([^,\t]+)/\1__\2/g' |   # convert last entry
  gzip -c >> "$outfile"

# compress
gzip -fn "$sintax_out"
