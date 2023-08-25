# Main configuration file

The configuration file lives in `config/config.yaml` of the analysis directory. It is structured as follows:

* **software**: pipeline-specific settings (version, paths, etc.)
* **input**: Paths to input file or file lists
* **workflows**: List of workflows that should be compared
* **primers**: List of primers to use (by barcode marker)
* **taxonomy_dbs**, *taxonomy_methods*: Taxonomy databases/assignment methods to use (see also [this documentation](taxonomy.md))
* **pipeline-specific settings**: Configuration options such as *cmp*, *filter*, *usearch*, *dada2*, *ITSx*. They define the default settings used by *all* workflows *unless* overridden specifically (see below).
  

## Examples

### USEARCH/VSEARCH-based pipeline (ITS)

```yaml
software:
  usearch:
    binary: usearch  # only needed with program: usearch

input:
  illumina:
    directories:
      - path/to/fastq_gz_files

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
    combinations: default
  trim_settings:
    min_overlap: 15
    max_error_rate: 0.25

# Settings for sequence processing/clustering
filter:
  min_length: 100

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
    ident_threshold: 0.97

# taxonomy assignment
taxonomy_dbs:
  ITS:
    unite:
      db: unite_all
      defined: order

taxonomy_methods:
  sintax_70:
    method: sintax_usearch
    confidence: 0.7

# ITS recognition (fungi)
ITSx:
  organism_groups: all
  e-value: 1.0e-5
  allow_single_domain: 1e-9,0
```

## Details on configuration sections

Disclaimer: This list is not complete, please refer to [this `config.yaml`](../config/config.yaml) for more information on some options.

### Input

Input files can be specified in various ways that best suit your needs, structured by sequencing technology (currently only *Illumina* possible). 

One thing is important to note: Paths are **always** relative to the *analysis directory* (not the *Snakecharmer* software directory, within which the pipeline is acutally run). If you want to be more clear about this, use absolute paths.

#### `directories`

The simplest is to just specify directories containing the sequence files (`fastq.gz`):

```yaml
input:
  illumina:
    directories:
      - path/to/fastq_directory
      - path/to/another_directory
    recursive: true  # searches subdirectories
    sample_pattern: illumina
```

Here, the setting `name_pattern` is of vital importance for extracting the sample names from the file names. You can either specify a known pattern (currently *illumina* for recent Illumina-style names and *sample_read* for something like `sample_name_R1.fastq.gz`). More complex cases need a regular expression pattern (more details in [this `config.yaml`](../config/config.yaml)).

#### `patterns`

A second option is to specify a pattern (or list of patterns) to specifically select certain files (e.g. a subset of a sequencing run, belonging to a specific dataset):

```yaml
input:
  illumina:
    patterns:
      - path/to/files/DatasetA*.fastq.gz
    sample_pattern: illumina
```

It is also possible to specify both `directories` and `patterns` (the selection may even overlap).

#### `sample_file`

If the above doesn't work for you, specify a sample file listing all sample names and paths. Example `samples.tsv`:

```
id  R1  R2
sample1 path/to/sample1_R1  path/to/sample1_R2
sample2 path/to/sample2_R1  path/to/sample2_R2
```

`config.yaml`:

```yaml
input:
  illumina:
    sample_file: samples.tsv
```

With single-end reads, just specify a two-column file without the `R2` column.

*Note*: [QIIME2 manifest files](https://docs.qiime2.org/2023.7/tutorials/importing/#fastq-manifest-formats) are accepted as well.

*Note 2*: The sample files are also produced internally if `directories` and/or `patterns` is specified. They are found here: `input/sample_config/<layout>/<run>/samples.tsv`

#### Sequencing runs

Normally, files are assumed to belong to the same sequencing run (called "run1" in file paths). If multiple runs of the same (or different) samples should be analyzed, there are two ways of specifying this (using `directories` in this example):

##### Listing each run

```yaml
input:
  illumina:
    run1:
      directories: [ path/to/run1 ]
      sample_pattern: illumina
    run1:
      directories: [ path/to/run2 ]
      sample_pattern: illumina
```

##### Using {run}

In the above example, we have a regular pattern (subdirectories residing next to each other), we could thus use the special `{run}` wildcard with the same result:

```yaml
input:
  illumina:
    directories: [ path/to/{run} ]
    sample_pattern: illumina
```

#### Run pooling

While it is usually a good idea to analyze each sequencing run separately (especially with DADA2), it can make sense to process data from different runs together in other cases (e.g. USEARCH-based clustering). The `pool_raw` setting will do exactly this. Samples are combined into a run pool (e.g. called `run1_run2_pool`), whereby reads from samples with the same name are combined. The default is not to pool, returing a separate results subdirectory for each run.

```yaml
input:
  illumina:
    directories: [ path/to/{run} ]
    sample_pattern: illumina

(...)

pool_raw: true
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
