# Amplicon data processing for diversity analysis using different workflows

This software makes use of the workflow management system [Snakemake](https://snakemake.github.io/) to build and integrate different amplicon workflows and allows for the flexible processing of data from different markers.

**Features**

- Denoising/clustering of raw amplicon sequencing data, taxonomic assignments and further sequence comparisons and marker-specific processing (such as ITS extraction)
- Comparison of different pipelines or variations of the same pipeline using different sets of parameters based on a flexible [`configuration system`](config/config.yaml). Downstream analyses and comparisons are easily possible due to a common file structure of the output directories.
- Simultaneous processing of multi-marker amplicons generated using different primer sets
- Multiple taxonomic assignment methods can be applied to each denoised dataset using marker-specific reference databases; currently implemented: [UNITE](https://unite.ut.ee) for Eukaryote ITS and [Midori](http://www.reference-midori.info) for mitochondrial markers (more may follow)

*Note:* To this date, a few workflows for paired-end Illumina data and taxonomic assignments for Eukaryote lineages been implemented, but the software may be extended further (see [below](#still-not-finished)).

**Non-features**

The software is mainly intended for the basic tasks of obtaining an OTU table, taxonomic assignments, sequence searches and basic validation of different workflows, even though the integrated pipelines often offer additional comprehensive downstream analyses. However, the output files are in commonly used formats (such as BIOM) and can serve as input for many analysis toolkits. There is also an [import script for R](#analyzing-in-r).

**Integrated pipelines**

- [USEARCH](https://www.drive5.com/usearch/manual)/[VSEARCH](https://github.com/torognes/vsearch)-based amplicon pipeline using UNOISE3 for obtaining ASVs (implemented here from scratch with Snakemake)
- [QIIME2](https://qiime2.org) (DADA2 denoising)
- [Amptk](https://github.com/nextgenusfs/amptk) (UNOISE3 and DADA2)
- ... (more to follow)

**Validation**

Validation is currently done using amplicon data from a fungal mock community ([details in `test` directory](test/README.md)) and the results of the different pipelines can be compared [using a script](#comparison-of-denoisingclustering-pipelines).

## Installing

The pipeline makes use of the [Conda package manager](https://conda.io), the installation is thus pretty straightforward ([see instructions here](INSTALL.md)).

## Configuring

The easiest is to copy the contents of the [config](config/) directory into the future analysis directory and then modify the files according to your needs. [`config.yaml`](config/config.yaml) contains all settings, while the taxonomic databases can be configured in [`taxonomy.yaml`](config/taxonomy.yaml). Please refer to the comments in these config files.

## Running

There are a few Snakemake target rules, which can be run independently (even though they also partly depend on each other). A complete list of commands is [found here](Commands.md).

### On a local computer

The following command runs the test pipeline (FASTQ files from sequencing fungal mock comunities in the [`test` directory](test/), specified with `-d test`) using 6 cores on a local computer. The  [target rules](Commands.md) to be run are `denoise`, `cmp`, `taxonomy` and `ITS`.

```sh
conda activate snakemake
snakemake -c6 --use-conda --conda-prefix ~/conda -d test denoise cmp taxonomy ITS
```

Note that the `~/conda` directory is used for the installation of all additional necessary software (`--conda-prefix` argument). This allows reusing the installed software across analyzes of different datasets.

### On a computer cluster

The command for processing a dataset on a HPCC may look like this (using SLURM with [this profile](https://github.com/Snakemake-Profiles/slurm#quickstart)):

```sh
outdir=~/path/to/analysis  # must contain a config directory
conda activate snakemake
snakemake -d $outdir \
    -j10 -c10 \
    --use-conda --conda-prefix ~/conda \
    --profile slurm \
    denoise cmp taxonomy ITS \
    --group-components sample=50 qc=100 \
    --rerun-incomplete \
    --default-resources runtime=180 mem_mb=500 \
    --rerun-triggers mtime
```

This command runs a maximum of 10 simultaneous jobs (`-j`), which may run on different nodes, each job with 20 cores (`-c`). Adding `--rerun-triggers mtime` at the end is currently recommended to avoid unnecessary re-running of QIIME/Amptk jobs.

#### Job grouping

For the USEARCH-based pipeline, it is recommended to use `--group-components` to limit the number of submitted jobs. Let's assume that 300 samples should be analyzed. By default, paired-end read merging, primer trimming and quality filtering is done separately for every of those samples on a single CPU core, meaning  that one job is submitted to the cluster for every sample. The same is true for the QC (FastQC) of the raw sequencing reads. However, `group-components` [allows processing](https://snakemake.readthedocs.io/en/v7.19.1/executing/grouping.html#job-grouping) multiple samples together; in the above example 50 of them (`sample=50`), resulting in only six sample processing jobs to be submitted instead of 300. Since the QC is done separately for forward and reverse reads, we group together 100 single read files (`qc=100`), which results in 6 QC jobs. After pre-processing, the sequences are combined for denoising, which can only run on a single node. However, if many pipelines/pipeline variants should be evaluated, component grouping can also be used (e.g. `--group-components denoise=4`).

Other pipelines (currently QIIME and Amptk) don't processing samples in parallel on different nodes. For these, the pre-processing steps belong to another group called `prepare`. Furthermore, the following groups exist: `otutab` (OTU table construction for USEARCH pipelines), `taxonomy` (taxonomy assignment), `ITS` (ITSx) and `cmp` (sequence comparisons). It is also possible to [assign single rules to custom groups](https://snakemake.readthedocs.io/en/v7.19.1/executing/grouping.html#job-grouping).

## Analyzing in R

The R source file [`R/read_amplicon.R`](R/read_amplicon.R) provides code for reading all data from a results directory. See also the small [example analysis](test/R_example/example.md).


## Comparison of denoising/clustering pipelines

There is a separate bash script `scripts/compare_results.sh`, which creates an Excel file comparing the number of reads assigned to each 98% cluster by each pipeline. A separate workbook is created for each sample. The script requires VSEARCH, as well as R with the following packages: `ggplot2`, `tidyverse`, `data.table` and `openxlsx`.


## Still not finished...

A list of possible next steps includes:

- Enable single-end analyses and other platforms than Illumina
- Integrate more pipelines / clustering methods
- Integrate more taxonomy databases
- Testing deployment on different systems
- Improve configuration of job resources (memory, CPUs)
- The USEARCH pipeline may be moved into an extra repository to be used independently. This repo focuses on the comparison of a variety of different pipelines as well as handling complex cases (e.g. several primer combinations).
- ...
