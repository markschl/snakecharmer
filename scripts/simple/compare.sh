#!/usr/bin/env bash

# This script compares the test pipeline results with
# simple example pipelines published by the pipeline authors
# in order to validate our code

set -xeuo pipefail

threads=12

# this relies on the 'simple' env
# conda env create -f scripts/simple_env.yaml

#conda activate simple

# prepare input files "manually"
gz=test/gz
out=test/processing/simple_input
rm -rf $out
mkdir -p $out
zcat $gz/mock1_R1.fastq.gz $gz/mock1_more/mock1_R1.fastq.gz > $out/mock1_R1.fastq
zcat $gz/mock1_R2.fastq.gz $gz/mock1_more/mock1_R2.fastq.gz > $out/mock1_R2.fastq
zcat $gz/mock2_R1.fastq.gz > $out/mock2_R1.fastq
zcat $gz/mock2_R2.fastq.gz > $out/mock2_R2.fastq

fprimers=$out/fprimers.fa
rprimers=$out/rprimers.fa
printf ">ITS3_KYO2\nGGGATGAAGAACGYAGYRAA\n>ITS3_KYO2\nTCGATGAAGAMCGYWGCVAD\n" > $fprimers
printf ">ITS4\nGCATATCAATAAGCGGAGGATT\n>ITS4\nGCATATTAWTCAGCGGAGGATT\n" > $rprimers

scripts/simple/vsearch.sh test file:$fprimers file:$rprimers $threads $out/*_R1.fastq
scripts/simple/usearch.sh test file:$fprimers file:$rprimers $threads $out/*_R1.fastq

scripts/compare_results.sh test

# VSEARCH pipeline presumably should yield same results
vres=test/results/unoise/data
ures=test/results/unoise_usearch/data
usearch_res=test/results/unoise_usearch_simple/pipeline_usearch_unoise3_simple/ITS__ITS3-KYO2...ITS4/paired
vsearch_res=test/results/unoise_vsearch_simple/pipeline_usearch_unoise3_simple/ITS__ITS3-KYO2...ITS4/paired
# case-insensitive sequence comparisons because of masking (our pipeline converts masked sequences to uppercase)
if ! cmp -s <(tr '[:upper:]' '[:lower:]' < $vres/denoised.fasta) <(tr '[:upper:]' '[:lower:]' < $vsearch_res/denoised.fasta); then
  echo "ASVs from simple vs. regular VSEARCH pipelines differ ($vres/denoised.fasta and $vsearch_res/denoised.fasta)" >&2
  exit 1
fi
if ! cmp -s $vres/denoised_otutab.txt.gz $vsearch_res/denoised_otutab.txt.gz; then
  echo "ASVs from simple vs. regular VSEARCH pipelines differ ($vres/denoised_otutab.txt.gz and $vsearch_res/denoised_otutab.txt.gz)" >&2
fi
# for USEARCH, direct comparisons are more difficult. This has to do with the de-replication,
# which our pipeline does separately per sample before collecting the unique sequences
# from all samples and de-replicating again. The standard USEARCH pipeline does it for all 
# samples together at once. The output is sorted by size, but in case of ties (same size),
# the order will vary between the two strategies. The output of U/VSEARCH clustering algorithms
# depends on the order of the input sequences, therefore the results could be different.
# assuming 'seqtool' is installed, we can compare sorted OTUs.
# For the OTU counts -> have a look at the Excel file in the 'cmp' directory
if ! cmp -s <(st . --to-tsv seq $ures/denoised.fasta | sort) <(st . --to-tsv seq $usearch_res/denoised.fasta | sort); then
  echo "ASVs from simple vs. regular USEARCH pipelines differ ($ures/denoised.fasta and $usearch_res/denoised.fasta)" >&2
  exit 1
fi

# for comparing biom:
# meld <(python -m json.tool $vres/denoised.biom) <(python -m json.tool $vsearch_res/denoised.biom)
# meld <(python -m json.tool $ures/denoised.biom) <(python -m json.tool $usearch_res/denoised.biom)
