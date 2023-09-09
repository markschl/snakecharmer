#!/bin/env bash

# this script contains the example VSEARCH pipeline
# https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline/c4859786f05bba35d8c306de4a3d64fea40d9dbf
# (adapted for UNOISE3, using the pipeline config)

set -euo pipefail

source scripts/parse_yaml.sh


if [ $# -lt 5 ]; then
    echo "usage: $0 <outdir> <f_primer> <r_primer> <threads> <fastq_files>..." 1>&2
    exit 1
fi

VSEARCH=vsearch

outdir="$1" && shift
f_primer="$1" && shift
r_primer="$1" && shift
threads="$1" && shift

out="$outdir"/workdir/amptk_simple
# obtain settings from config file
eval $(parse_yaml "$outdir/config/config.yaml")

rm -rf "$out"
mkdir -p "$out"

fq_dir="$out/fq"
mkdir -p "$fq_dir"
cp "$@" "$fq_dir"

set -x

########### start ##################################

conda activate $amptk

# calculate absolute number of primer mismatches
pmismatch=$(bc <<< "(${#f_primer} + ${#r_primer})/2 * $primers_trim_settings_max_error_rate" | printf %.0f -)

# merge/trim
amptk illumina -i "$fq_dir" -o "$out" \
   -f "$f_primer" -r "$r_primer" \
    -f {params.f_primer_seq} -r {params.r_primer_seq} \
    --min_len $filter_min_length \
    --trim_len 10000000 `# high enough to never be longer`  \
    --cpus $threads \
    --cleanup \
    --require_primer=on \
    --rescue_forward=off \
    --primer_mismatch $pmismatch \
    --usearch $(which usearch)

# cluster


################################################################

# finish up
# Copy to results dir in order to enable a comparison with the other pipelines

out="$outdir"/results/unoise_vsearch_simple/pipeline_usearch_unoise3_simple/ITS__ITS3-KYO2...ITS4/paired
mkdir -p $out
cp otus.fasta $out/clusters.fasta
gzip -cn otutab.txt > $out/otutab.txt.gz
