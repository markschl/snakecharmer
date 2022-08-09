#!/usr/bin/env bash

set -euo pipefail



if [ $# -lt 3 ]; then
  echo "Usage: $0 <infile> <refdb> <outfile> [sintax_args]" >&2
  exit 1
fi


infile="$1" && shift
refdb_compr="$1" && shift
outfile="$1" && shift

out_name=$(basename $outfile)
outdir=$(dirname $outfile)
sintax_out=$outdir/sintax/"${out_name%.*}"
refdb="${outfile%.gz}"_ref.fasta

mkdir -p $outdir/sintax

# first, uncompress the database
zstd -dqf "$refdb_compr" -o "$refdb"

usearch -sintax "$infile" -db "$refdb" \
  -tabbedout "$sintax_out" \
  -quiet \
  -log /dev/stdout \
  "$@" 2>&1

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
gzip -f "$sintax_out"
