#!/usr/bin/env bash

set -xeuo pipefail


if [ $# -lt 2 ]; then
    echo "usage: $0 <url_prefix> <outfile>" >&2
    echo "  url: file URL prefix (without '.fasta.gz' or '.taxon.gz'" >&2
    echo "  outfile: *.fasta.zst" >&2
    exit 1
fi

url_prefix="$1" && shift
outfile="$1" && shift

if [[ ! "$url_prefix" =~ 'QIIME' ]]; then
  echo "QIIME not found in Midori file URL, please specify QIIME-formatted version" >&2
  exit 1
fi

# make temporary directory
tmp_dir="${outfile%.*.*}"_tmp
mkdir -p "$tmp_dir"

# obtain files
echo "Downloading..." >&2
fa="$tmp_dir/midori.fasta"
tax="$tmp_dir/midori.taxon"
wget -q --no-check-certificate -O "$fa.gz" "$url_prefix.fasta.gz"
wget -q --no-check-certificate -O "$tax.gz" "$url_prefix.taxon.gz"
gzip -df "$fa.gz" "$tax.gz"

echo "Converting..." >&2

if [[ "$url_prefix" =~ /QIIME/ ]]; then
  # Filtered version without sp. / aff. / ..., only sequences with species label present
  # compose FASTA file with taxonomy
  st set -ul "$tax" -d '{l:2}' "$fa" |
    zstd -c > $outfile

elif [[ "$url_prefix" =~ /QIIME_sp/ ]]; then
  # Non-filtered version including sp. and others defined only at higher ranks
  # -> we have to convert undefined names to empty strings.
  # The patterns follow information from https://doi.org/10.1002/edn3.303 and comparisons with
  # filtered files from Midori. Still, it was not possible to obtain the exact
  # same result, differences are however very small (see validation below).
  sp='sp\.|sp\.?_\d+|aff\.|cf\.|gen\.|nom\.|nud\.|_nomen .*|_nom\. .*|nr\.|BOLD|unidentified|undescribed|uncultured|complex|unclassified|taxon \w+'

  # compose FASTA file with taxonomy and convert unknown ranks to empty string
  st set -ul "$tax" -d '{l:2}' "$fa" |  # add taxonomy to description
    st replace -rd ";g__genus_.*" ';g__;s__' |  # remove undefined genus + species if g__genus is found (some hybrids retained by Midori are removed, unfortunately)
    st replace -rd ";s__([^ ]+ +)*?($sp)( [^;]+)?$" ';s__' |  # remove species matching one of the above patterns
    st replace -rd ";f__family_.+?;g__;" ";f__;g__;" |  # remove undefined family, etc.
    st replace -rd ";o__order_.+?;f__;" ";o__;f__;" |
    st replace -rd ";c__class_.+?;o__;" ";c__;o__;" |
    st replace -rd ";p__phylum_.+?;c__;" ";p__;c__;" |
    st replace -d ' ' '_' |
    zstd -c > $outfile

  # validate by comparing the filtered version by Midori
  filt=$(sed 's/_sp//gi' <<< "$url_prefix")
  tax2="$tmp_dir/midori_filtered.taxon"
  if wget -q --no-check-certificate -O "$tax2.gz" "$filt.taxon.gz"; then
    gzip -df "$tax2.gz"
    sort "$tax2" | tr ' ' '_' > "$tax2.sorted"
    zstd -dcq "$outfile" | 
      st find -fr --desc 's__[^;]+' --to-tsv id,desc |  # exclude taxa defined at species level (after above replacements)
      sort \
      > "$tax.species.sorted"
    diff "$tax2.sorted" "$tax.species.sorted" > "$tmp_dir/comparison.txt" || true
    n0=$(< "$tax2" wc -l)
    n=$(< "$tax" wc -l)
    nfilt=$(< "$tax.species.sorted" wc -l)
    echo "$n sequences found, $nfilt have a valid species label (according to this script), while $n0 are valid species according to Midori." >&2
    n0=$(grep -e '^<' $tmp_dir/comparison.txt | wc -l)
    n1=$(grep -e '^>' $tmp_dir/comparison.txt | wc -l)
    echo "$n0 sequences are only in the filtered Midori sequence set, but not recognized as valid species by this script." >&2
    echo "$n1 sequences recognized by this script as valid species, but not by Midori" >&2
  else
    echo "Filtered Midori database could not be downloaded, no comparison possible" >&2
  fi
fi

# clean up
rm -rf "$tmp_dir"
