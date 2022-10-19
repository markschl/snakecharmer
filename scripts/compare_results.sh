#!/usr/bin/env bash

# TODO: this script may be incorporated in the pipeline

set -euo pipefail


if [ $# -eq 0 ]; then
  echo "Usage: $0 <project_dir>" >&2
  echo "This script compares the denoising results of different pipelines." >&2
  echo "Output: Excel files in cmp directory" >&2
  exit 1
fi

project_dir=$1

# settings

# clustering identities
id=0.98
# min. identity for mapping of ASVs against common clusters
map_id=0.97
# min. identity for self comparison of clusters
self_map_id=0.9
# number of CPUs to use
cpu=4

mkdir -p $project_dir/cmp
echo -n > $project_dir/cmp/cmp.txt
echo -n > $project_dir/cmp/self_cmp.txt

# collect ASVs/OTUS by primer combination
ls $project_dir/results/*/*/*/*/denoised.fasta | 
  sed -E "s|$project_dir/([^/]+/){3,3}([^/]+).*|\2|g" |
  sort -u |
  while read primers; do
    echo "$primers"
    # cluster all the OTUs/ASVs together at defined threshold
    cluster_file=$project_dir/cmp/clusters.$primers.fasta
    cat $project_dir/results/*/*/$primers/*/denoised.fasta |
      vsearch -cluster_fast - \
        -id $id \
        -centroids $cluster_file \
        -relabel seq  \
        -threads $cpu \
        -maxaccepts 64 -maxrejects 64

    # map every OTU set against these clusters
    for f in $project_dir/results/*/*/$primers/*/denoised.fasta; do
      strategy=$(basename $(dirname $f))
      name=$(basename $(dirname $(dirname $(dirname $(dirname $f)))))
      vsearch -usearch_global $f -db $cluster_file \
        -id $map_id -maxhits 1 -quiet \
        -userout - -userfields query+target+id \
        -threads $cpu \
        -maxaccepts 64 -maxrejects 64 |
        awk -F$'\t' -v OFS=$'\t'  \
          -v name=$name -v primers=$primers -v strategy=$strategy \
          '$0=name FS primers FS strategy FS $0' >> $project_dir/cmp/cmp.txt
    done

    # also, map clusters against themselves to obtain the closest relatives
    vsearch -usearch_global $cluster_file -db $cluster_file \
      -id $self_map_id -self -quiet \
      -userout - -userfields query+target+id \
      -threads $cpu \
      -maxaccepts 64 -maxrejects 64 |
      awk -F$'\t' -v OFS=$'\t' -v primers=$primers \
        '$0=primers FS $0' >> $project_dir/cmp/self_cmp.txt
  done


scripts/compare_results.R $project_dir

rm $project_dir/cmp/cmp.txt $project_dir/cmp/self_cmp.txt
