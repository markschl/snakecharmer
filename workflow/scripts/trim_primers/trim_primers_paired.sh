#!/usr/bin/env bash

set -xeuo pipefail

if [ $# -lt 4 ]; then
    echo "usage: $0 <input_file> <fwd_primers> <rev_primers_rev> <outdir> [cutadapt options]" 1>&2
    echo "Note: the reverse primers need to be reverse complemented..." 1>&2
    exit 1
fi

input_file="$1" && shift
forward="$1" && shift
reverse_rev="$1" && shift
outdir="$1" && shift

sample=$(basename ${input_file%.fastq.zst})
# note: multiqc does not yet parse cutadapt.json, so we use simple log files

# remove old files in case there are some
rm -f "$outdir/"*...*.fastq
rm -f "$outdir/"*...*.fastq.zst

# recognize forward and reverse primers
zstd -dcq "$input_file" |
  cutadapt - \
    -g "file:$forward" \
    --suffix ' fwd={name}' \
    `# --json $outdir/"$sample"_fwd.cutadapt.json` \
    "$@" 2> $outdir/"$sample"_fwd.log |
  cutadapt - \
    -a "file:$reverse_rev" \
    --suffix ' rev={name}' \
    `# --json $outdir/"$sample"_rev.cutadapt.json` \
    "$@" 2> $outdir/"$sample"_rev.log |
  st split --fq -o "$outdir/{a:fwd}...{a:rev}.fastq"

zstd -q --rm "$outdir/"*...*.fastq

# add sample name to cutadapt logs: a hack needed for multiqc to recognize the sample
sed -i -E "s/(Command line parameters[^$]+$)/\1 $sample.fastq.gz/g" $outdir/"$sample"_fwd.log
sed -i -E "s/(Command line parameters[^$]+$)/\1 $sample.fastq.gz/g" $outdir/"$sample"_rev.log

# parse logfiles to obtain total numbers
log=$outdir/"$sample"_rev.log
stats="$outdir/$sample"_stats.txt
n=$(grep 'Total reads processed' $outdir/"$sample"_fwd.log | sed -E 's/[^0-9]+([0-9,]+)/\1/g' | tr -d ',')
printf "$n" > $stats
for dir in fwd rev; do
  n_trimmed=$(grep 'Reads with adapters' $outdir/"$sample"_"$dir".log | sed -E 's/[^0-9]+([0-9,]+) *\(.+/\1/g' | tr -d ',')
  printf "\t$n_trimmed" >> $stats
done
