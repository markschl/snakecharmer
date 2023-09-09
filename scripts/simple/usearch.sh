#!/usr/bin/env bash

# This script follows https://www.drive5.com/usearch/manual/ex_miseq_its.html
# (https://www.drive5.com/usearch/manual/ex_miseq_its.bash)
# as closely as possible to reproduce the workflow and 
# match the results against those of our pipeline
# (the approach is the same)

set -euo pipefail

source scripts/simple/parse_yaml.sh


if [ $# -lt 5 ]; then
    echo "usage: $0 <outdir> <f_primer> <r_primer> <threads> <fastq_files>..." 1>&2
    exit 1
fi

usearch=usearch

outdir="$1" && shift
f_primer="$1" && shift
r_primer="$1" && shift
threads="$1" && shift

out="$outdir"/workdir/usearch_simple
# obtain settings from config file
eval $(parse_yaml "$outdir/config/config.yaml")

rm -rf "$out"
mkdir -p "$out"

set -x

########### start ##################################

# Merge all R1/R2 file pairs
# Add sample name to read labels (-relabel @ option)
# Pool samples together into raw.fq
$usearch -fastq_mergepairs "$@" -relabel @ -fastqout "$out"/raw.fq \
   -fastq_pctid $usearch_merge_overlap_ident \
   -fastq_maxdiffs $usearch_merge_max_diffs \
   -threads $threads
# fastq_pctid and fastq_maxdiffs added based on pipeline config

########### extra code ##########################

# Primer trimming is needed with this data
# We trim forward/reverse primers separately to ensure that only trimmed
# reads are retained in case of primer mixes (linked adapters would not be possible?)
cutadapt "$out/raw.fq" \
    -g "$f_primer" \
    --error-rate $primers_trim_settings_max_error_rate \
    --overlap $primers_trim_settings_min_overlap \
    --cores $threads \
    --discard-untrimmed \
    -o "$out/trimmed_fwd.fq"

cutadapt "$out/trimmed_fwd.fq" \
    -a "$r_primer" \
    --error-rate $primers_trim_settings_max_error_rate \
    --overlap $primers_trim_settings_min_overlap \
    --minimum-length $usearch_filter_min_length \
    --cores $threads \
    --discard-untrimmed \
    -o "$out/trimmed.fq"

# now we can cd in to the directory to retain the code as unchanged as possible
cd "$out"

##################################################

# Quality filter
$usearch -fastq_filter trimmed.fq -fastq_maxee_rate $usearch_filter_max_error_rate \
  -fastaout filtered.fa -relabel Filt \
  -threads $threads
# using maxee_rate instead of maxee (presumably better with variable length amplicons)

# Find unique read sequences and abundances
$usearch -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa

# Run UPARSE algorithm to make 97% OTUs
$usearch -cluster_otus uniques.fa -otus otus.fa -relabel Otu \
    -threads $threads

# Run UNOISE algorithm to get denoised sequences (ZOTUs)
$usearch -unoise3 uniques.fa -zotus zotus.fa \
  -minsize $usearch_unoise3_min_size \
  -threads $threads
# (minsize added based on pipeline config)

# Downstream analysis of OTU sequences & OTU table
# Can do this for both OTUs and ZOTUs, here do
# just OTUs to keep it simple.
##################################################

# Make OTU table
$usearch -otutab trimmed.fq -otus zotus.fa -otutabout otutab_raw.txt \
  -id $(bc <<< "scale=4; $usearch_otutab_ident_threshold / 100") \
  -maxaccepts $usearch_otutab_maxaccepts \
  -maxrejects $usearch_otutab_maxrejects \
  -biomout otutab_raw.biom \
  -threads $threads
# (using zotus.fa instead of otus.fa; maxaccepts, maxrejects added based on pipeline config)

# # Normalize to 5k reads / sample
# $usearch -otutab_norm otutab_raw.txt -sample_size 5000 -output otutab.txt

# # Alpha diversity
# $usearch -alpha_div otutab.txt -output alpha.txt

# # Make OTU tree
# $usearch -calc_distmx otus.fa -tabbedout distmx.txt
# $usearch -cluster_aggd distmx.txt -treeout otus.tree

# # Beta diversity
# mkdir beta/
# $usearch -beta_div otutab.txt -tree otus.tree -filename_prefix beta/

# # Rarefaction
# $usearch -alpha_div_rare otutab.txt -output alpha_rare.txt

# # Predict taxonomy
# $usearch -sintax otus.fa -db ../data/rdp_its_v2.fa -strand both \
#   -tabbedout sintax.txt -sintax_cutoff 0.8

# # Taxonomy summary reports
# $usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank g -output genus_summary.txt
# $usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank p -output phylum_summary.txt


################################################################

# finish up

out=../../results/unoise_usearch_simple/workflow_usearch_unoise3_simple/run1_run2_pool_paired/ITS__ITS3-KYO2...ITS4
mkdir -p $out
cp zotus.fa $out/denoised.fasta
gzip -cn otutab_raw.txt > $out/denoised_otutab.txt.gz
cp otutab_raw.biom $out/denoised.biom
