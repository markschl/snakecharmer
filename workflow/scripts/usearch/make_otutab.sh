#!/usr/bin/env bash

set -xeuo pipefail

# FIXME: this script should be converted to a regular Snakemake script (Python or Bash),
# it has far too many arguments
if [ $# -lt 8 ]; then
    echo "usage: $0 <program> <threads> <uniques> <otus> <otutab_out> <map_out> <bam_out> <notmatched_out>  [usearch options]" 1>&2
    exit 1
fi

program="$1" && shift
threads="$1" && shift
uniques_compr="$1" && shift
otus="$1" && shift
tab="$1" && shift
map="$1" && shift
bam="$1" && shift
notmatched="$1" && shift

sam=${bam%.*}.sam
notmatched=${notmatched%.*}

if [[ "$program" == "vsearch" ]]; then
    zstd -dcq "$uniques_compr" |
    vsearch -usearch_global - -db "$otus" \
        -otutabout "${tab%.gz}" \
        -userout "${map%.gz}" \
        -userfields query+target+id \
        -maxhits 1 \
        -strand plus \
        -sizein \
        -samout "$sam" \
        -notmatched "$notmatched" \
        -threads $threads \
        "$@"

    # SAM -> BAM
    rm -f "$otus".fai "$bam".bai
    samtools view -T "$otus" -b $sam |
        samtools sort -@ $threads > "$bam"
    rm $sam "$otus".fai
    samtools index "$bam"

elif [[ "$program" == "usearch" ]]; then
    uniques=$(mktemp ${uniques_compr%.fasta.zst}.XXXXXX.fasta)
    zstd -dqf "$uniques_compr" -o "$uniques"
    # TODO: -mapout does not include the identity
    usearch -otutab "$uniques" -otus "$otus" \
        -otutabout "${tab%.gz}" \
        -userout "${map%.gz}" \
        -userfields query+target+id \
        -maxhits 1 \
        -notmatched "$notmatched" \
        -threads $threads \
        "$@" 1>&2
    rm "$uniques"

    # SAM produced by USEARCH is sometimes not parsable
    # -> create "stub" BAM file just to make sure it is there
    echo -n > "$bam"
else
    echo "unknown program: {params.par[program]}"
    exit 1
fi
# compress
gzip -n "${tab%.gz}" "${map%.gz}"
zstd -qf "$notmatched"
