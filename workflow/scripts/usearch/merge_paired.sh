#!/usr/bin/env bash

set -xeuo pipefail

# FIXME: this script should be converted to a regular Snakemake script (Python or Bash),
# it has far too many arguments
if [ $# -lt 6 ]; then
    echo "usage: $0 <forward>.fastq.gz <reverse>.fastq.gz <out_prefix> <program> <usearch_bin> <min_ident_pct> [usearch options]" 1>&2
    exit 1
fi

forward="$1" && shift
reverse="$1" && shift
out_prefix="$1" && shift
program="$1" && shift
usearch_bin="$1" && shift
min_ident_pct="$1" && shift

fwd="$out_prefix"_R1.fastq
rev="$out_prefix"_R2.fastq
merged="$out_prefix".fastq
nm_fwd="$out_prefix"_notmerged_R1.fastq
nm_rev="$out_prefix"_notmerged_R2.fastq

# first, extract the input files
gzip -dc "$forward" > "$fwd"
gzip -dc "$reverse" > "$rev"

# perform the merging and parse the output (different for programs)

if [[ "$program" == "vsearch" ]]; then
    max_diffpct=$( bc <<< "100 - $min_ident_pct" )
    merge_out=$( \
        vsearch \
            --fastq_mergepairs "$fwd" \
            --reverse "$rev" \
            --fastqout "$merged" \
            --fastqout_notmerged_fwd "$nm_fwd" \
            --fastqout_notmerged_rev "$nm_rev" \
            --fastq_allowmergestagger \
            --fastq_maxdiffpct "$max_diffpct" \
            "$@" 2>&1 \
        )
    
elif [[ "$program" == "usearch" ]]; then
    merge_out=$( \
        "$usearch_bin" \
            -fastq_mergepairs "$fwd" \
            --reverse "$rev" \
            -fastqout "$merged" \
            -fastqout_notmerged_fwd "$nm_fwd" \
            -fastqout_notmerged_rev "$nm_rev" \
            -fastq_pctid "$min_ident_pct" \
            "$@" 2>&1 \
        )
else
  echo "Read merging: program argument needs to be either 'usearch' or 'vsearch'" >&2
  exit 1
fi

# output stats file based on parsing of intercepted usearch output
# (fortunately, the output of USEARCH and VSEARCH is similar enough for a common parsing routine)
n=$(grep " Pairs" <<< "$merge_out" | sed -E 's/ *([0-9]+) *Pairs.*/\1/g')
n_merged=$(grep " Merged" <<< "$merge_out" | sed -E 's/ *([0-9]+) *Merged.*/\1/g')
printf "$n\t$n_merged" > "$out_prefix"_stats.txt

# clean up
rm "$fwd" "$rev"

# finally, compress the files
zstd -q --rm "$merged"
zstd -q --rm "$nm_fwd"
zstd -q --rm "$nm_rev"
