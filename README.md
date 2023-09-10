# Amplicon data processing workflows for diversity analysis

This software makes use of the workflow management system [Snakemake](https://snakemake.github.io/) to build and integrate different amplicon workflows and allows for the flexible processing of data from different markers.

**Features**

- *Clustering/denoising* of raw amplicon sequencing data, *taxonomic assignments* and further sequence comparisons and *marker-specific processing* (such as ITS extraction)
- *Comparison* of different pipelines or variations of the same pipeline using *different sets of parameters*. Workflow output is presented in a common file structure.
- Simultaneous processing of *multi-marker amplicons* generated using different primer sets
- *Multiple taxonomic assignment methods* can be applied to each dataset using different marker-specific reference databases

This is a work in progress and will be extended further (see [below](#further-steps)).

**Non-features**

The software does not assist with comprehensive statistical analyses, even though the integrated pipelines may offer them. Since the output files are in commonly used formats (such as BIOM), they can still serve as input for many analysis toolkits. There is also a dedicated [import script for R](#analyzing-in-r).

**Integrated pipelines**

- [uvsnake](https://github.com/markschl/uvsnake), an [USEARCH](https://www.drive5.com/usearch/manual)/[VSEARCH](https://github.com/torognes/vsearch)-based paired-end amplicon pipeline (UPARSE and UNOISE3)
- [QIIME2](https://qiime2.org) with DADA2 denoising offering paired-end and single-end analyses
- [Amptk](https://github.com/nextgenusfs/amptk) (UNOISE3 and DADA2, paired-end only)
- ... (more to follow)

**Validation**

Validation is done using amplicon data from a fungal mock community ([details in `test` directory](test/README.md)) and a basic comparison of the different workflows can be done [using a dedicated script](#comparison-of-pipelines).

No snakes 🐍 were harmed in the process of creating this software

## Installing

The software makes use of the [Conda package manager](https://conda.io), the installation is thus pretty straightforward ([see instructions here](INSTALL.md)).

## Configuring

The easiest is to copy the contents of the [config](config/) directory into a new analysis directory and rename `config.template.yaml` to `config.yaml` and `taxonomy.template.yaml` to `taxonomy.yaml`, and then modify the files according to your needs. The two configuration files:

* **`config.yaml`**: Main configuration file containing input and workflow definitions. [See here](docs/config.md) for a (incomplete) description. The available options are are further documented [in `config/config.template.yaml`](config/config.template.yaml).
* **`taxonomy.template.yaml`**: Defines all available taxonomcic databases. [See here for details](docs/taxonomy.md), as well as [`config/taxonomy.template.yaml`](config/taxonomy.template.yaml) for examples.

## Running

The `snakecharmer` script is used as follows:

```
./snakecharmer <outdir> <rule1> <rule2>...
```

Before running the first time, the command `conda activate snakemake` is necessary if [Conda was used](INSTALL.md).

The target rules (or "commands") define, which output should be generated. Only one or several rules can be specified. Some depend on output of other rules. For instance, the `taxonomy` rule requires the clustering/denoising to happen before (`cluster` command). A complete list of commands is [found here](docs/rules.md).

The most important are:

* `quality`: Reports sequencing read statistics in `results/<workflow>/validation`, which may help in deciding on setting the workflow parameters.
* `cluster`: Does the clustering/denoising for all workflows defined in `config.yaml`.
* `taxonomy`: Applies all taxonomy assignment methods defined in `config.yaml` to the output of all clustering workflows.
* `clean`: Removes working directories that are not strictly needed (retaining the `results` dir)

### Example

The following command processes a test dataset (fungal mock comunities in the [`test` directory](test/)) using 6 cores on a local computer. The  [target rules](docs/rules.md) to be run are `cluster`, `cmp`, `taxonomy` and `ITS`.

```sh
conda activate snakemake
./snakecharmer test cluster cmp taxonomy  --cores 6
```

On a computer cluster, the command may look different ([see documentation here](https://snakemake.readthedocs.io/en/latest/executing/cluster.html)). The `snakecharmer` script is just a simple wrapper for Snakemake, but otherwise accepts all arguments that Snakemake does. Example for running with SLURM:

```sh
outdir=~/path/to/analysis
./snakecharmer --cores 20 --jobs 20 --slurm $outdir cluster cmp taxonomy
```

### Output

After running, a few additional directories will have appeared next to `config`. The most important one is the `results` directory, which roughly has the following structure ([more details here](docs/output.md)):

```
📦<my_analysis>/
 ├─ 📂 config/
 │  ├─ 🗋 config.yaml
 │  └─ 🗋 taxonomy.yaml
 │  (...)
 ├─ 📂 results/
 │  ├─ 📂 <workflow name>/
 │  │  ├─ 📂 data/
 │  │  │  ├─ 🗋 clusters.fasta
 │  │  │  ├─ 🗋 otutab.txt.gz
 │  │  │  ├─ 🗋 otutab.biom
 │  │  │  ├─ 🗋 otutab.hdf5.biom
 │  │  │  ├─ 📂 taxonomy/
 │  │  │  │  ├─ 🗋 <database>-<method>-<name>..txt.gz
 │  │  │  │  ├─ 🗋 <database>-<method>-<name>.biom.gz
 │  │  │  │  ├─ 🗋 <database>-<method>-<name>.hdf5.biom.gz
 │  │  │  │  │  (...)
 │  │  │  ├─ 📂 cmp/
 │  │  │  │  ├─ 🗋 <my_seq_comparison>.txt
 │  │  │  │  ├─ 🗋 <my_seq_comparison>_notmatched.fasta.gz
 │  │  │  │  ├─ 🗋 <my_seq_comparison>_clusters_notmatched.fasta.gz
 │  │  │  │  │  (...)
 │  │  │  ├─ 📂 [ITSx/]
 │  │  │  │  ├─ 🗋 out.positions.txt
 │  │  │  │  └─ (...)
 (...)
```

Whether the output are ASVs/ESVs or OTUs from a fixed threshold clustering (not yet implemented), the resulting FASTA file is always called `clusters.fasta`. The sample/OTU count matrix is returned both in the traditional tabular format (`otutab.txt.gz`) and [BIOM](https://biom-format.org/documentation/biom_format.html). The taxonomic annotations are named by taxonomy database and assignment method (multiple combinations possible) and returned in a QIIME-style tabular format as well as the BIOM format. Furthermore, there can be results of sequence comparisons (`cmp`) or marker-specific data such as ITSx results.

With the simplest scenario (one run/layout and one primer combination), the relevant results directory is `<my_analysis>/results/<workflow_name>/data`. With multi-workflow/marker setups, the `data` directory will not be present, and the individual workflow results are placed in the nested directories ([more here](docs/output.md)).

## Analyzing in R

The R source file [`R/read_amplicon.R`](R/read_amplicon.R) provides code for reading all data from a results directory. See also the small [example analysis](test/R_example/example.md).

## Comparison of denoising/clustering pipelines

There is a separate bash script `scripts/compare_results.sh`, which creates an Excel file comparing the number of reads assigned to 98% clusters of the already clustered sequences by each pipeline. A separate workbook is created for each sample. The script requires VSEARCH, as well as R with the following packages: `ggplot2`, `tidyverse`, `data.table` and `openxlsx`.

## Further steps...

A list of possible next steps includes:

- Integrate more pipelines / clustering methods and taxonomy databases
- Integrate other platforms than Illumina and allow simultaneous analysis of multi-platform data
- Allow for run result merging
- Offer more ways of comparing and validating pipelines and generally improve user experience
- Testing deployment on different systems
- Improve configuration of job resources (memory, CPUs)
- ...

