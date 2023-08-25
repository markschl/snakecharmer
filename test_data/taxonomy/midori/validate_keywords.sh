#!/usr/bin/env bash

# Code for validating filter keywords
wget -O - https://www.reference-midori.info/forceDownload.php?fName=download/Databases/GenBank256/QIIME/longest/MIDORI2_LONGEST_NUC_GB256_CO1_QIIME.taxon.gz | zcat | sort -u > filtered
wget -O - https://www.reference-midori.info/forceDownload.php?fName=download/Databases/GenBank256/QIIME_sp/longest/MIDORI2_LONGEST_NUC_SP_GB256_CO1_QIIME.taxon.gz | zcat | sort -u > all

# entries only in the non-filtered dataset
comm -23 all filtered > all_only
# entries only in the filtered dataset (there are few, must be some inconsistency)
comm -13 all filtered > filtered_only
# exclude all sequences matching the given species pattern
# (trying to be as close as possible to the Midori filtered version)
# (also removing g__genus_ entries (but without g__genus_[Name]), 
# since the pipeline script also sets species from these to undefined)
# The pattern in the python script is on multiple lines and was condensed here, with \s replaced by spaces
pat="[_ ]((cf\.|aff\.|sp\.|environment|undescribed|uncultured|complex|unclassified|nom\.|_nomen .*|_nom\. .*|nud\.|unidentif\.|indet\.|gen\.|nr\.|taxon \w+)[_ \d]|(sp|cf)[\._][A-Z0-9]|sp\.$)"
grep -Ev "$pat" all | grep -Ev "g__genus_[^\[]" > all_filtered

# only in Midori filtered version, missing from ours
# (meaning we may have filtered too restrictively)
# We find mostly species with g__genus_, from which a majority are hybrids, which may
# make sense to exclude anyway. <10 names are real species that for some reason
# don't have a genus set. These are unfortunately set to undefined.
comm -13 all_filtered filtered  # 87

# only in our filtered version but missing from the Midori version
# (Midori might have been too restrictive there)
comm -23 all_filtered filtered  # 21

# or, use two versions downloaded using this pipeline
# (using https://github.com/markschl/seqtool)
f1=filtered.fasta.zst
f2=all.fasta.zst
zstdcat $f1 | st . --to-tsv id,desc | sort -u > filtered
zstdcat $f2 | st . --to-tsv id,desc |sort -u > all
comm -23 all filtered > all_only
comm -13 all filtered > filtered_only
