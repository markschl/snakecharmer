#!/usr/bin/env bash

set -xeuo pipefail


if [ $# -lt 3 ]; then
    echo "usage: $0 <input_file> <sample_name> <out_prefix> [vsearch-filter-args]" 1>&2
    exit 1
fi

input_file="$1" && shift
sample_name="$1" && shift
out_prefix="$1" && shift

#### filter ####

out=$( {  
    zstd -dcq "$input_file" |  # decompress file
      st set -i "$sample_name.{num}" --fq |  # set sequence number and delete description
      st del -d --fq |  # delete descriptions
      vsearch -fastq_filter - \
        -fastaout - \
        -fastaout_discarded "$out_prefix"_discarded.fasta \
        "$@" \
      > "$out_prefix"_filter.fasta; 
  } 2>&1 )

# output stats file based on parsing of intercepted output of VSEARCH command
grep 'sequences kept' <<< "$out" | 
  sed -E 's/ *([0-9]+)[^,]+, *([0-9]+).*/\1 \2/g' | 
  tr ' ' '\t' > "$out_prefix"_stats.txt

# forward output
cat <<< "$out" 1>&2


#### de-replicate ####

# "good" uniques only used for clustering

vsearch -derep_fulllength "$out_prefix"_filter.fasta \
    -sizeout \
    -output "$out_prefix"_good_uniques.fasta

# all uniques: used for constructing the OTU table (note -sizein argument)

cat "$out_prefix"_good_uniques.fasta "$out_prefix"_discarded.fasta |
  vsearch -derep_fulllength - \
      -sizein -sizeout \
      -output - |
    zstd -cq > "$out_prefix"_all_uniques.fasta.zst

# compress
zstd -fq --rm "$out_prefix"_filter.fasta
zstd -fq --rm "$out_prefix"_good_uniques.fasta
# TODO: always produced by usearch?
zstd -fq --rm "$out_prefix"_discarded.fasta
