#!/usr/bin/env bash

set -xeuo pipefail

if [ $# -lt 4 ]; then
    echo "usage: $0 <input_file> <fwd_primers> <rev_primers_rev> <outdir> <min_length> [cutadapt options]" 1>&2
    echo "Note: the reverse primers need to be reverse complemented..." 1>&2
    exit 1
fi

input_file="$1" && shift
forward="$1" && shift
reverse_rev="$1" && shift
outdir="$1" && shift
min_length="$1" && shift

sample=$(basename ${input_file%.fastq.zst})
short_file="$outdir/too_short.fastq"
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
    --minimum-length $min_length \
    --too-short-output "$short_file" \
    `# --json $outdir/"$sample"_rev.cutadapt.json` \
    "$@" 2> $outdir/"$sample"_rev.log |
  st split --fq -o "$outdir/{a:fwd}...{a:rev}.fastq"

if [ -e "$short_file" ]; then
  zstd -qf --rm "$short_file"
fi

# add sample name to cutadapt logs: a hack needed for multiqc to recognize the sample
sed -i -E "s/(Command line parameters[^$]+$)/\1 $sample.fastq.gz/g" $outdir/"$sample"_fwd.log
sed -i -E "s/(Command line parameters[^$]+$)/\1 $sample.fastq.gz/g" $outdir/"$sample"_rev.log

# statistics
if grep -q "No reads processed!" $outdir/"$sample"_fwd.log; then
  n=0
  n_trimmed_f=0
  n_trimmed_r=0
  n_long=0
else
  # parse logfiles to obtain total numbers
  extract_num() {
    sed -E 's/[^0-9]+([0-9,]+).*/\1/g' | tr -d ','
  }
  n=$(grep 'Total reads processed' "$outdir/$sample"_fwd.log | extract_num)
  n_trimmed_f=$(grep 'Reads with adapters' "$outdir/$sample"_fwd.log | extract_num)
  n_trimmed_r=$(grep 'Reads with adapters' "$outdir/$sample"_rev.log | extract_num)
  n_long=$(grep 'Reads written' "$outdir/$sample"_rev.log | extract_num)
  # TODO (when turning into Python script): this reports only the reverse sequences that are long enough,
  # not very intuitive
fi
printf "$n\t$n_trimmed_f\t$n_trimmed_r\t$n_long" > "$outdir/$sample"_stats.txt
