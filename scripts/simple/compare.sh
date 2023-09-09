#!/usr/bin/bash

# This script compares the test pipeline results with
# simple example pipelines published by the pipeline authors
# in order to validate our code

set -euo pipefail

# In order for results to be reproducible, we use only one core.
# This is especially important for USEARCH, where the order of 
# merged reads changes with multithreading (apparently not with VSEARCH).
threads=1

# this relies on the 'uvsearch' env
# conda env create -f scripts/simple/uvsearch_env.yaml

source $(conda info --base)/etc/profile.d/conda.sh
conda activate uvsearch

# prepare input files "manually"
gz=test/gz
out=test/processing/simple_input
rm -rf $out
mkdir -p $out
zcat $gz/run1/mock1_R1.fastq.gz $gz/run2/mock1_R1.fastq.gz > $out/mock1_R1.fastq
zcat $gz/run1/mock1_R2.fastq.gz $gz/run2/mock1_R2.fastq.gz > $out/mock1_R2.fastq
zcat $gz/run1/mock2_R1.fastq.gz > $out/mock2_R1.fastq
zcat $gz/run1/mock2_R2.fastq.gz > $out/mock2_R2.fastq

fprimers=$out/fprimers.fa
rprimers=$out/rprimers.fa
printf ">ITS3_KYO2\nGGGATGAAGAACGYAGYRAA\n>ITS3_KYO2\nTCGATGAAGAMCGYWGCVAD\n" > $fprimers
printf ">ITS4\nGCATATCAATAAGCGGAGGATT\n>ITS4\nGCATATTAWTCAGCGGAGGATT\n" > $rprimers

# Run VSEARCH analysis
scripts/simple/vsearch.sh test file:$fprimers file:$rprimers $threads $out/*_R1.fastq
# Run USEARCH analysis
scripts/simple/usearch.sh test file:$fprimers file:$rprimers $threads $out/*_R1.fastq

scripts/compare_results.sh test

# VSEARCH pipeline presumably should yield same results
vres=test/results/unoise/data
ures=test/results/unoise_usearch/data
usearch_res=test/results/unoise_usearch_simple/workflow_usearch_unoise3_simple/run1_run2_pool_paired/ITS__ITS3-KYO2...ITS4
vsearch_res=test/results/unoise_vsearch_simple/workflow_usearch_unoise3_simple/run1_run2_pool_paired/ITS__ITS3-KYO2...ITS4
# case-insensitive sequence comparisons because of masking (our pipeline converts masked sequences to uppercase)
if ! cmp -s <(tr '[:upper:]' '[:lower:]' < $vres/clusters.fasta) <(tr '[:upper:]' '[:lower:]' < $vsearch_res/clusters.fasta); then
  echo "ASVs from simple vs. regular VSEARCH pipelines differ ($vres/clusters.fasta and $vsearch_res/clusters.fasta)" >&2
  exit 1
fi
if ! cmp -s $vres/otutab.txt.gz $vsearch_res/otutab.txt.gz; then
  echo "ASVs from simple vs. regular VSEARCH pipelines differ ($vres/otutab.txt.gz and $vsearch_res/otutab.txt.gz)" >&2
fi
# for USEARCH, direct comparisons are more difficult. This has to do with the de-replication,
# which our pipeline does separately per sample before collecting the unique sequences
# from all samples and de-replicating again. The standard USEARCH pipeline does it for all 
# samples together at once. The output is sorted by size, but in case of ties (same size),
# the order will vary between the two strategies. The output of U/VSEARCH clustering algorithms
# depends on the order of the input sequences, therefore the results could be different.
# assuming 'seqtool' is installed, we can compare sorted OTUs.
# For the OTU counts -> have a look at the Excel file in the 'cmp' directory
if ! cmp -s <(st . --to-tsv seq $ures/clusters.fasta | sort) <(st . --to-tsv seq $usearch_res/clusters.fasta | sort); then
  echo "ASVs from simple vs. regular USEARCH pipelines differ ($ures/clusters.fasta and $usearch_res/clusters.fasta)" >&2
  exit 1
fi

# for comparing biom with Meld:
# meld <(python -m json.tool $vres/otutab.biom) <(python -m json.tool $vsearch_res/otutab.biom)
# meld <(python -m json.tool $ures/otutab.biom) <(python -m json.tool $usearch_res/otutab.biom)
