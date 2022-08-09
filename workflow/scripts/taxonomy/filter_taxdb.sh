#!/usr/bin/env bash

# Downloads a reference database from a specified location and creates a compressed
# FASTA file containing sequence IDs and QIIME-formatted taxonomy as description line.
# The taxonomy will have the whole lineage defined (no ranks missing). If names
# are not known, then they should be empty.

set -euo pipefail


if [ $# -lt 3 ]; then
    echo "usage: $0 <infile> <filter> <outfile>" >&2
    echo "  - infile: zstd-compressed FASTA file in QIIME taxonomy format" >&2
    echo "  - filter options: domain, kingdom, class, phylum, order, family or species" >&2
    echo "  - outfile: filtered output file, zstd-compressed" >&2
    exit 1
fi

infile="$1" && shift
filter="$1" && shift
outfile="$1" && shift

if [[ "$filter" =~ ^(species|genus|family|order|class|phylum|kingdom)$ ]]; then

  echo "Filter taxonomy: keeping sequences defined at $filter rank in $outfile" >&2

  # create filter regex expression
  # e.g. f__Family or s__Some_species should be present and not empty
  char=${filter::1}  # character code
  expr="$char"__"[^\s$;]+"  # filter expression with preceding semicolon

  # regex filtering with seqtool
  dropped="${outfile%.*.*}"_dropped.fasta
  zstd -dc "$infile" |
    st find -fr --desc "$expr" --dropped "$dropped" |
    zstd -qc > "$outfile"
  
  # check if numbers match (if seqtool did it correctly)
  # TODO: comment
  orig=$(zstd -dc "$infile" | st count)
  n_out=$(zstd -dc "$outfile" | st count)
  n_removed=$(st count "$dropped")
  if [[ $(( $n_out + $n_removed )) -ne $orig ]]; then
    echo "Sequence numbers before/after filtering don't match: $n_out + $n_removed != $orig" >&2
    exit 2
  fi
  # clean up
  # TODO: don't check dropped seqs at all
  rm "$dropped"

else
  echo "Unknown taxonomy filter: $filter" >&2
  exit 1

fi
