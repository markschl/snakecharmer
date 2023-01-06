#!/bin/env bash

# this script contains the example VSEARCH pipeline
# https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline/c4859786f05bba35d8c306de4a3d64fea40d9dbf
# (adapted for UNOISE3, using the pipeline config)

set -euo pipefail

source scripts/simple/parse_yaml.sh


if [ $# -lt 5 ]; then
    echo "usage: $0 <outdir> <f_primer> <r_primer> <threads> <fastq_files>..." 1>&2
    exit 1
fi

VSEARCH=vsearch

outdir="$1" && shift
f_primer="$1" && shift
r_primer="$1" && shift
THREADS="$1" && shift

out="$outdir"/processing/vsearch_simple
# obtain settings from config file
eval $(parse_yaml "$outdir/config/config.yaml")
set -x

rm -rf "$out"
mkdir -p "$out"


########### start ##################################

$VSEARCH --fastq_chars $(ls -1 "$@" | head -1)

# Process samples                                                               

for f in "$@"; do

    r=$(sed -e "s/_R1/_R2/" <<< "$f")
    s=$(basename "$f" | cut -d_ -f1)

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================
    echo

    $VSEARCH --fastq_mergepairs $f \
        --threads $THREADS \
        --reverse $r \
        --fastq_maxdiffs $usearch_merge_max_diffs \
        --fastq_maxdiffpct $(bc <<< "100 - $usearch_merge_overlap_ident") \
        --fastq_allowmergestagger \
        --fastqout "$out"/$s.merged.fastq \
        --fastq_eeout
    # Added fastq_allowmergestagger, fastq_maxdiffs, fastq_maxdiffpct based on pipeline config

    # Commands to demultiplex and remove tags and primers
    # using e.g. cutadapt may be added here.

    ################ extra code #########################################

    # Primer trimming is needed with this data
    # We trim forward/reverse primers separately to ensure that only trimmed
    # reads are retained in case of primer mixes (linked adapters would not be possible?)
    cutadapt "$out"/$s.merged.fastq \
        -g "$f_primer" \
        --error-rate $primers_trim_settings_max_error_rate \
        --overlap $primers_trim_settings_min_overlap \
        --cores $THREADS \
        --discard-untrimmed \
        -o "$out"/$s.trimmed_fwd.fastq

    cutadapt "$out"/$s.trimmed_fwd.fastq \
        -a "$r_primer" \
        --error-rate $primers_trim_settings_max_error_rate \
        --overlap $primers_trim_settings_min_overlap \
        --minimum-length $filter_min_length \
        --cores $THREADS \
        --discard-untrimmed \
        -o "$out"/$s.trimmed.fastq

    # Addition: Prepare a FASTA version with the sample name to allow constructing
    # an OTU table from unfiltered reads
    $VSEARCH --fastx_filter "$out"/$s.trimmed.fastq \
      --relabel $s. \
      --fastaout "$out"/$s.trimmed.fasta

    ######################################################################

    echo
    echo Calculate quality statistics
    echo

    $VSEARCH --fastq_eestats "$out"/$s.trimmed.fastq \
        --output "$out"/$s.stats
    echo
    echo Quality filtering
    echo

    $VSEARCH --fastq_filter "$out"/$s.trimmed.fastq \
        --fastq_maxee_rate $usearch_filter_max_error_rate \
        --fastq_minlen $filter_min_length \
        --fastaout "$out"/$s.filtered.fasta \
        --fasta_width 0
    ## these options from the example workflow are not used:
        # --fastq_maxlen 275 \
        # --fastq_maxns 0 \
    ## additionally using fastq_maxee_rate

    echo
    echo Dereplicate at sample level and relabel with sample.n
    echo

    $VSEARCH --derep_fulllength "$out"/$s.filtered.fasta \
        --strand plus \
        --output "$out"/$s.derep.fasta \
        --sizeout \
        --relabel $s. \
        --fasta_width 0

done

# Now we can cd in to the output directory to keep the code
# as similar as possible
cd "$out"

# At this point there should be one fasta file for each sample
# It should be quality filtered and dereplicated.

echo
echo ====================================
echo Processing all samples together
echo ====================================
echo
echo Merge all samples

cat *.derep.fasta > all.fasta

echo
echo Sum of unique sequences in each sample: $(cat all.fasta | grep -c "^>")
echo
echo Dereplicate across samples
echo

$VSEARCH --derep_fulllength all.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.derep.uc \
    --output derep.fasta

# echo
# echo Unique sequences: $(grep -c "^>" derep.fasta)
# echo
# echo Cluster sequences using VSEARCH
# echo

# $VSEARCH --cluster_size derep.fasta \
#     --threads $THREADS \
#     --id $CLUSTERID \
#     --strand plus \
#     --sizein \
#     --sizeout \
#     --fasta_width 0 \
#     --centroids centroids.fasta

#echo
#echo Cluster with Swarm using d=1 and fastidious mode
#echo
#
#$SWARM derep.fasta \
#    --threads $THREADS \
#    --differences 1 \
#    --fastidious \
#    --seeds centroids.fasta \
#    --usearch-abundance \
#    --output /dev/null

# Instead of the above clustering commands, implement UNOISE3
# done as described in https://github.com/torognes/vsearch/pull/283
$VSEARCH --cluster_unoise derep.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --minsize $usearch_unoise_min_size \
    --centroids zotus_chim.fasta

$VSEARCH --sortbysize zotus_chim.fasta \
    --output zotus_sorted.fasta
    
$VSEARCH --uchime3_denovo zotus_sorted.fasta \
    --nonchimeras centroids.fasta

echo
echo Clusters: $(grep -c "^>" centroids.fasta)
echo
echo Sort and remove singletons
echo

$VSEARCH --sortbysize centroids.fasta \
    --threads $THREADS \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --output sorted.fasta
# options from example workflow not used:
    # --minsize 2

echo
echo Non-singleton clusters: $(grep -c "^>" sorted.fasta)
echo 
echo De novo chimera detection
echo

# Using uchime3_denovo instead of uchime_denovo
# (difference is abskew=16 instead of 2, https://github.com/torognes/vsearch/pull/283)
# (we omit uchime_ref)
$VSEARCH --uchime3_denovo sorted.fasta \
    --sizein \
    --nonchimeras otus.fasta \
    --relabel Zotu
# options from example workflow not used:
#  --fasta_width 0 \
#  --sizeout \
# --qmask none

# echo
# echo Unique sequences after de novo chimera detection: $(grep -c "^>" denovo.nonchimeras.fasta)
# echo
# echo Reference chimera detection
# echo

# $VSEARCH --uchime_ref denovo.nonchimeras.fasta \
#     --threads $THREADS \
#     --db $REF \
#     --sizein \
#     --sizeout \
#     --fasta_width 0 \
#     --qmask none \
#     --dbmask none \
#     --nonchimeras nonchimeras.fasta

# echo
# echo Unique sequences after reference-based chimera detection: $(grep -c "^>" nonchimeras.fasta)
# echo
# echo Relabel OTUs
# echo

# $VSEARCH --fastx_filter nonchimeras.fasta \
#     --threads $THREADS \
#     --sizein \
#     --sizeout \
#     --fasta_width 0 \
#     --relabel OTU_ \
#     --fastaout otus.fasta


echo
echo Number of OTUs: $(grep -c "^>" otus.fasta)
echo
echo Map sequences to OTUs by searching
echo

$VSEARCH --usearch_global all.fasta \
    --threads $THREADS \
    --db otus.fasta \
    --id $usearch_otutab_ident_threshold \
    --strand plus \
    --sizein \
    --sizeout \
    --maxaccepts $usearch_otutab_maxaccepts \
    --maxrejects $usearch_otutab_maxrejects \
    --fasta_width 0 \
    --qmask none \
    --dbmask none \
    --otutabout otutab.txt
# added maxaccepts and maxrejects based on the pipeline config 

# We also construct an OTU table based on all (unfiltered) sequences,
# as recommended by the USEARCH author and currently done by default in our pipeline
cat *.trimmed.fasta >  trimmed.fasta
$VSEARCH --usearch_global trimmed.fasta \
    --threads $THREADS \
    --db otus.fasta \
    --id $usearch_otutab_ident_threshold \
    --strand plus \
    --sizein \
    --sizeout \
    --maxaccepts $usearch_otutab_maxaccepts \
    --maxrejects $usearch_otutab_maxrejects \
    --qmask none \
    --dbmask none \
    --otutabout otutab.txt \
    --biomout otutab.biom


# echo
# echo Sort OTU table numerically
# echo

# sort -k1.5n otutab.txt > otutab.sorted.txt

echo
echo Done

################################################################

# finish up
# Copy to results dir in order to enable a comparison with the other pipelines

out=../../results/unoise_vsearch_simple/pipeline_usearch_unoise3_simple/ITS__ITS3-KYO2...ITS4/paired
mkdir -p $out
cp otus.fasta $out/denoised.fasta
gzip -cn otutab.txt > $out/denoised_otutab.txt.gz
cp otutab.biom $out/denoised.biom
