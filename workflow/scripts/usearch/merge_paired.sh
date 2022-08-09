#!/usr/bin/env bash

set -xeuo pipefail


if [ $# -lt 3 ]; then
    echo "usage: $0 <forward>.fastq.gz <reverse>.fastq.gz <out_prefix> [usearch options]" 1>&2
    exit 1
fi

forward="$1" && shift
reverse="$1" && shift
out_prefix="$1" && shift

fwd="$out_prefix"_R1.fastq
rev="$out_prefix"_R2.fastq
merged="$out_prefix".fastq
nm_fwd="$out_prefix"_notmerged_R1.fastq
nm_rev="$out_prefix"_notmerged_R2.fastq

# first, extract the input files
gzip -dc "$forward" > "$fwd"
gzip -dc "$reverse" > "$rev"

# then perform the merging, intercepting the output in a variable
# TODO: make usearch/vsearch configurable
out=$(usearch \
    -fastq_mergepairs "$fwd" \
    --reverse "$rev" \
    -fastqout "$merged" \
    -fastqout_notmerged_fwd "$nm_fwd" \
    -fastqout_notmerged_rev "$nm_rev" \
    "$@" 2>&1)

# clean up
rm "$fwd" "$rev"

# output stats file based on parsing of intercepted usearch output
n=$(grep Pairs <<< "$out" | sed -E 's/ *([0-9]+) *Pairs.*/\1/g')
n_merged=$(grep Merged <<< "$out" | sed -E 's/ *([0-9]+) *Merged.*/\1/g')
printf "$n\t$n_merged" > "$out_prefix"_stats.txt

# forward usearch output
cat <<< "$out" 1>&2

# finally, compress the files
zstd -q --rm "$merged"
zstd -q --rm "$nm_fwd"
zstd -q --rm "$nm_rev"
