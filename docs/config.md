# Main configuration file

The configuration file lives in `config/config.yaml` of the analysis directory. It is structured as follows:

* **software**: pipeline-specific settings (version, paths, etc.)
* **input**: Paths to input file or file lists
* **workflows**: List of workflows that should be compared
* **primers**: List of primers to use (by barcode marker)
* **taxonomy_dbs**, *taxonomy_methods*: Taxonomy databases/assignment methods to use (see also [this documentation](taxonomy.md))
* **pipeline-specific settings**: Configuration options such as *cmp*, *filter*, *usearch*, *dada2*, *ITSx*. They define the default settings used by *all* workflows *unless* overridden specifically (see below).
  

## Examples

### USEARCH/VSEARCH-based pipeline

This example uses UNOISE3 clustering of an ITS2 amplicon and assigns taxonomy using SINTAX.

```yaml
software:
  usearch:
    binary: usearch  # only needed with program: usearch
  uvsnake:
    snakemake_env: snakemake
    repo:
      tag: v0.1

input:
  illumina:
    run1: path/to/run1_sample_tab.tsv
    run2: path/to/run2_sample_tab.tsv

workflows:
  # UNOISE3 pipeline using VSEARCH for denoising
  unoise:
    cluster: usearch_unoise3
    taxonomy: default
  # UNOISE3 (VSEARCH) pipeline using a different filtering threshold
  unoise_strict:
    cluster: usearch_unoise3
    taxonomy: default
    settings:
      usearch:
        filter: { max_error_rate: 0.001 }
  # UNOISE3 pipeline using USEARCH for read merging
  unoise_umerge:
    cluster: usearch_unoise3
    taxonomy: default
    settings:
      usearch:
        merge: { program: usearch }

primers:
  ITS:
    forward: 
      - ITS3-KYO2: GGGATGAAGAACGYAGYRAA
    reverse:
      - ITS4: AATCCTCCGCTTATTGATATGC
  trim_settings:
    min_overlap: 15
    max_error_rate: 0.25
    min_length: 100

# Settings for sequence processing/clustering
usearch:
  # paired end merging
  merge:
    overlap_ident: 75  # percent
    max_diffs: 1000
    expected_length: 400

  filter:
    max_error_rate: 0.002  #  0.8 errors per 400 bp

  unoise:
    min_size: 8

  otutab:
    maxaccepts: 64
    maxrejects: 64
    ident_threshold: 97  # percent

# taxonomy assignment
taxonomy:
  dbs:
    ITS:
      unite:
        db: unite_all
        defined: order

  methods:
    sintax_70:
      method: sintax_usearch
      confidence: 0.8

# ITS recognition (fungi)
ITSx:
  organism_groups: all
  e-value: 1.0e-5
  allow_single_domain: 1e-9,0
```

## Details on configuration sections

Disclaimer: This list is not complete, please refer to [`config/config.template.yaml`](../config/config.template.yaml) for more information on some options.

### Input

Input files are be specified as sample lists with read file paths. Each list has two or three columns with the sample name, forward and (optionally) reverse paths. [QIIME2 manifest files](https://docs.qiime2.org/2023.7/tutorials/importing/#fastq-manifest-formats) are accepted as well.

A simple way of creating a sample table given one or several directories is the [`make_sample_tab`](https://github.com/markschl/ngs-sample-tab) script. Given this directory structure with Illumina sequencing reads (names were already simplified):

```
raw_reads/
├── run1/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
├── run2/
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
│   ├── sample3_R1.fastq.gz
│   ├── sample3_R2.fastq.gz
```

The following will create two sample files. 

```sh
make_sample_tab --recursive --directory raw_reads --format simple
# or: make_sample_tab -rd raw_reads -f simple
```

This yields `samples_run1_paired.tsv`:

```
id	R1	R2
sample1	raw_reads/run1/sample1_R1.fastq.gz	raw_reads/run1/sample1_R2.fastq.gz
sample2	raw_reads/run1/sample2_R1.fastq.gz	raw_reads/run1/sample2_R2.fastq.gz
```

...and `samples_run2_paired.tsv`:

```
id	R1	R2
sample2	raw_reads/run2/sample2_R1.fastq.gz	raw_reads/run2/sample2_R2.fastq.gz
sample3	raw_reads/run2/sample3_R1.fastq.gz	raw_reads/run2/sample1_R2.fastq.gz
```

The corresponding entry in `config/config.yaml`:

```yaml
input:
  illumina:
    run1: samples_run1_paired.tsv
    run2: samples_run2_paired.tsv
```

#### Run pooling

While it is usually a good idea to analyze each sequencing run separately (especially with DADA2), it can make sense to process data from different runs together in other cases (e.g. USEARCH-based clustering). The `pool_raw` setting will do exactly this. Samples are combined into a run pool (e.g. called `run1_run2_pool`), whereby reads from samples with the same name are combined. The default is not to pool, returing a separate results subdirectory for each run. In the above example, `sample2` would be pooled due to overlapping names.

```yaml
input:
  illumina:
    run1: samples_run1_paired.tsv
    run2: samples_run2_paired.tsv

(...)

pool_raw: true
```

Alternatively, samples could be made unique with `--unique`:

```sh
make_sample_tab --unique --pattern 'pool=raw_reads/run*/*.fastq.gz' --format simple
# short: make_sample_tab -up 'pool=raw_reads/run*/*.fastq.gz' -f simple 
```

This yields:

This yields `samples_pool_paired.tsv`:

```
id	R1	R2
sample1	raw_reads/run1/sample1_R1.fastq.gz	raw_reads/run1/sample1_R2.fastq.gz
sample2_1	raw_reads/run1/sample2_R1.fastq.gz	raw_reads/run1/sample2_R2.fastq.gz
sample2_2	raw_reads/run2/sample2_R1.fastq.gz	raw_reads/run2/sample2_R2.fastq.gz
sample3	raw_reads/run2/sample3_R1.fastq.gz	raw_reads/run2/sample3_R2.fastq.gz
```

## Primers

Forward and reverse amplicon primers have to be listed for them to be trimmed before further processing. The trimming settings are added in the same section and used by all pipelines. Currently, the UVSnake and QIIME2 pipelines make use of all settings, while we limit mismatches to 2 for Amptk ([explanation in code](../workflow/scripts/amptk_trim_paired.py)).

Primers are grouped per marker. The marker name needs to be present in the [`taxonomy`/`dbs`](taxonomy.md) configuration as well. Only databases from the given marker will be used for classification. Furthermore, some rules are marker-specific. For example, 'ITSx' only runs
on ITS amplicons (make sure to always use `ITS` as marker name).

*Example:*

```yaml
primers:
  marker:
    forward: 
      - fwd_primer: SEQUENCE
    reverse:
      - rev_primer: SEQUENCE
  trim_settings:
    min_overlap: 15
    max_error_rate: 0.25
    min_length: 100
```

### Primer mixes

Each primer can further be a comma delimited list of oligos, which were mixed together. Such mixes can be advantageous over simple degenerate primers to limit the number of variants present in the mix (see e.g. [Tedersoo 2015](https://doi.org/10.3897/mycokeys.10.4852)).

The oligos should all still start at the same nucleotide position to allow global trimming.

```yaml
primers:
  marker:
    forward: 
      - fwd_primer:
        SEQUENCE1,
        SEQUENCE2,
        SEQUENCE3
    reverse:
      - rev_primer:
        SEQUENCE1,
        SEQUENCE2
  (...)
```

**Important note:**: From the integrated pipelines some (QIIME2, Amptk) cannot deal with primer mixes out of the box, only [UVSnake](https://github.com/markschl/uvsnake) can. For the others, we currently use the **consensus sequence** of all primers in the mix.
The consensus threshold is 50%, meaning that at each position, 50% of the bases are represented in the consensus. This may lead to some rare ambiguous variants being removed, while adding some other unneeded wobbles. A higher threshold consensus threshold could lead to too many unspecific hits. With a sufficiently high `max_error_rate`, primers should still be trimmed well.
More options on this may be added in the future.

### Anchoring

The oligos may include adapter and/or linker sequences, but their position is not
required to be at the exact sequence start or end. If this is necessary,
you can anchor the sequences using: `^SEQUENCE`. UVSnake and QIIME2 will currently
use this setting.

### Combinations

Multiple primers can be specified in each direction. All combinations will be searched (per marker):


```yaml
primers:
  marker1:
    forward: 
      - fwd_primer1: FORWARD1
      - fwd_primer2: FORWARD2
    reverse:
      - rev_primer1: FORWARD1
      - rev_primer2: FORWARD2
  marker2:
    forward:
      - fwd_primer3: FORWARD3
    reverse:
      - rev_primer3: FORWARD3
  (...)
```

After trimming, the clustered sequences, OTU tables, etc. are found in the following directory:
`results/<workflow>/workflow_<method>/<run>/marker__<fwd primer>...<rev primer>/`
(see also [this documentation](output.md)).

Combinations can also be limited as follows:


```yaml
primers:
  marker1:
    forward: 
      - fwd_primer1: FORWARD1
      - fwd_primer2: FORWARD2
    reverse:
      - rev_primer1: FORWARD1
      - rev_primer2: FORWARD2
    combinations:
      - fwd_primer1...rev_primer1
      - fwd_primer2...rev_primer1
      - fwd_primer2...rev_primer2
  (...)
```

## Workflows

An unlimited number of workflows with parameters variations can be listed. Each workflow needs an unique name, and the output will be found in `results/<workflow name>`.

The general (minimal) structure is as follows:

```yaml
workflows:
  workflow1:
    cluster: cluster_method_a
    taxonomy: default
  workflow2:
    cluster: cluster_method_b
    taxonomy: default
  ...

cluster_method_a_settings:
  some_threshold: 0.5

cluster_method_b_settings:
  some_setting: some_value
```

This runs `workflow1` and `workflow2`, which use the clustering/denoising pipelines "a" and "b" with default settings listed below. Extra detail: some advanced settings are actually not required to be listed in `config.yaml`. All these settings and their default values can currently only be found in this [JSON Schema file](../workflow/config.schema.yaml).

### Cluster methods

The following clustering methods are available:

* *uvsnake* pipeline:
  * `uvsnake_unoise3`
  * `uvsnake_uparse`
* *QIIME2*:
  * `qiime2_dada2`
* *Amptk*:
  * `amptk_unoise3`
  * `amptk_uparse`
  * `amptk_dada2`

### Parameter variations

In order to modify the parameters for a specific workflow, use the `settings` section:

```yaml
workflows:
  workflow1_default:
    cluster: cluster_method_a
    taxonomy: default
  workflow1_30:
    cluster: cluster_method_a
    taxonomy: default
    settings:
      cluster_method_a_settings:
        some_threshold: 0.3
  workflow1_70:
    cluster: cluster_method_a
    taxonomy: default
    settings:
      cluster_method_a_settings:
        some_threshold: 0.7
  workflow1_pooled:
    settings:
      pool_raw: true
  ...
```

Within `settings`, just specify the settings you want to modify using the same hierarchy. You can see whether the default threshold was correctly changed by looking at `results/<workflow>/config.yaml`. For example, `results/workflow1_30/config.yaml` should contain the the correct setting for `some_threshold`:

```yaml
name: workflow1_30
cluster: cluster_method_a
taxonomy:
  (...)
settings:
cluster_method_a_settings:
  some_threshold: 0.3
(...)
```
