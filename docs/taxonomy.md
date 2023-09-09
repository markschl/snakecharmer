# Taxonomic databases

There are two ways of specifying databases in `config/taxonomy.yaml`:

* Sequence databases from various local and online sources.
  The exact format varies depending on the source (`format` setting). The data
  is automatically converted to FASTA files with QIIME-style taxonomy annotations
  internally, and can be further filtered and converted to the right format
  used by the classifier. Depending on the algorithm, the training may be
  resource-intensive and has the potential to use all memory of your local computer.
* Download pre-trained/pre-formatted databases.
  These can be used directly with the corresponding taxonomic assignment
  program.
  On the downside, they cannot be filtered or converted to other 
  formats and are thus tied to one or few specific assignment algorithms.

## Configuration

Here a simple/typical `taxonomy.yaml` with ITS (UNITE), COI (Midori) and 16S (SILVA) databases (more details in [this `config/taxonomy.template.yaml`](../config/taxonomy.template.yaml)):

```yaml
taxonomy_db_sources:
  unite_all:
    format: unite_otus
    doi: 10.15156/BIO/2483918
  Midori_COI:
    format: midori
    version: 256
    marker: lrRNA
    kind: longest
    remove_num: true
  SILVA_V4_QZA:
    format: qiime_qza
    sequences: https://data.qiime2.org/2023.5/common/silva-138-99-seqs-515-806.qza
    taxonomy: https://data.qiime2.org/2023.5/common/silva-138-99-tax-515-806.qza
    # some extra settings for recognizing/cleaning ambiguous names
    ambig:
      - undefined
      - unidentified
      - uncultured
    ambig_sp:
      - regex: "sp\\.?$"
```

In `config.yaml`, the databases are referenced by their name.

```yaml
taxonomy_dbs:
  ITS:
    unite:
      db: unite_all
    unite_order:
      db: unite_all
      defined: order  # drops sequences not annotated at the order level
  COI:
    midori:
        db: Midori_COI
  16S:
    silva:
        db: SILVA_V4_QZA
```

Taxonomy assignment is done using all databases available for a given barcode/marker. In this example, the ITS marker has `unite` and `unite_order`, allowing to compare the effect of dropping "poorly" annotated sequences from the reference database.

Markers are defined in the *primers* setting in `config.yaml`.


## Database formats

### General formats

A few general FASTA database formats are supported (local or remote):

* **qiime**: flat FASTA files with QIIME-style headers (also used internally in this
    pipeline). For example, UNITE [provides ITS databases](https://unite.ut.ee/repository.php)
    in this format (but see also *unite_otus* below).
* **qiime_qza**: QZA artifact archives used as database format by QIIME 2,
   e.g. obtained [from here](https://docs.qiime2.org/2023.5/data-resources)
   or prepared with [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt)
* **utax**: UTAX-formatted FASTA files, to be used with *sintax_usearch* 
   assignment method (see *taxonomy_methods* in config.yaml).
   Some files are provided on the
   [USEARCH home page](https://www.drive5.com/usearch/manual/sintax_downloads.html)

The `ambig` and `ambig_sp` setting further allow specifying keywords or patterns
recognizing and removing ambiguous/undefined names. The required settings vary
from source to source.

### Specific online sources

Converters are adapted to specific properties of these sources, where possible
ambiguous/undefined names are known.
The converters *unite_otus* and *midori* allow selecting the files by DOI or
other settings without knowing the exact URL.

* **unite_otus** for ITS (https://unite.ut.ee/repository.php),
  referenced by `DOI`
* **midori** for eukaryote mitochondrial markers 
   (https://www.reference-midori.info/download.php),
   referenced by `version`, `marker` and `kind`
* **gtdb** for 16S (https://gtdb.ecogenomic.org), referenced by URL(s)

*note*: it is still possible to download these databases manually and provide
   a file path.

## Pre-trained datasets

Pre-trained databases can be obtained from local or remote sources.
and may be gzip-compressed (usually .gz, .gzip).
Tarfiles or ZIP archives containing several files should be manually extracted
and the correct file supplied using path=...
Be careful not to use pre-trained databases from untrusted sources.


* **qiime_nb**: use with classifier from QIIME2, `qiime_sklearn` method.
   e.g. https://docs.qiime2.org/2023.5/data-resources
* **idtaxa**, use with `idtaxa`` assignment method;
   http://www2.decipher.codes/Downloads.html


## Assembling your own database

One way is to use [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt).

Here some commands from their [tutorial](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494),
downloading and processing data from SILVA (which is not directly possible
by this pipeline).

```sh
# obtain
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

# U -> T
qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza 
    --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

# low-quality
qiime rescript cull-seqs \
    --i-sequences silva-138.1-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

# length filter
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.1-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza 

# dereplicate
qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.1-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza
```

The following entry in `config/taxonomy.yaml` uses our new database:

```yaml
taxonomy_db_sources:
  silva_138_unique:
    format: qiime_qza
    sequences: path/to/silva-138.1-ssu-nr99-seqs-derep-uniq.qza
    taxonomy: path/to/silva-138.1-ssu-nr99-tax-derep-uniq.qza
```
