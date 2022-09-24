# Analysis of environmental diversity using amplicon sequencing

This pipeline system uses [Snakemake](https://snakemake.github.io/) to create and integrate several different amplicon pipelines. The results are returned in a common directory structure and common file formats and may serve as input for R-based analyses. The configuration should be easy, with a flexible system for specifying the input files and an automatic recognition of sample names and paired-end vs. single-end sequencing strategies.

Currently, only paired-end Illumina data can be analyzed, and the analyses are tailored for fungal amplicon (taxonomy assignments using UNITE).

Pipelines:

- [USEARCH](https://www.drive5.com/usearch/manual)/[VSEARCH](https://github.com/torognes/vsearch)-based amplicon pipeline using UNOISE3 for obtaining ASVs (implemented here from scratch with Snakemake)
- [QIIME2](https://qiime2.org) (DADA2)
- The [Amptk pipeline](https://github.com/nextgenusfs/amptk) (UNOISE3 and DADA2 denoising strategies)

The current rule graph:

![rule graph](rulegraph.png)


## Installing

The Conda package manager is the most important prerequisite. See here for [installation instructions in the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

Qiime is installed using Conda according to [these instructions](https://docs.qiime2.org/2022.8/install/native/#install-qiime-2-within-a-conda-environment). Make sure that the correct QIIME version is set in [`config/config.yaml`](config/config.yaml).

For Amptk, see [these instructions](https://amptk.readthedocs.io/en/latest/#install) (again using Conda). In addition, there is currently a [bug](https://github.com/nextgenusfs/amptk/issues/96) with UNOISE3 denoising. The code can be fixed like this:

```sh
conda activate amptk
sed -i "s/--derep_fulllength', filter_out, '--relabel', 'Read_', '--sizeout', '--output/--fastx_uniques', filter_out, '--relabel', 'Read_', '--sizeout', '--fastaout/g" "$CONDA_PREFIX/lib/python3.10/site-packages/amptk/unoise3.py"
```

All further software is automatically installed when running the pipeline (make sure to specify `--use-conda`).

An exception is *Usearch*. It can be [obtained here](https://www.drive5.com/usearch/download.html) and must be accessible as `usearch` in `$PATH`.

## Configuring

The easiest approach is to copy the contents of the [config](config/) directory into the future analysis directory and then modify the files according to your needs. `config.yaml` contains all settings, while the taxonomic databases can be configured in `taxonomy.yaml`.

Important options:
- `input`: both Bash glob patterns and/or input directories (optional recursive search) can be specified. Matching files will be linked in an `input` directory, grouped into single/paired-end sequencing strategies. If the same filename is encountered twice, they are currently merged together. A list of input samples is placed in `results/samples.yaml`
- `pipelines`: An unlimited list of independent pipelines with different clustering (denoising) methods can be specified. In addition, every entry can modify the default settings, which are also in `config.yaml` below the pipelines definition. To check, which settings were used by a pipeline see `results/<pipeline>/config.yaml` after starting the execution (or `snakemake -j1 config`).

## Running

This command runs the test pipeline on FASTQ files in the `test` directory (see [test](test/)) on a local computer, using `~/conda` for the installation of all additional necessary software:

```sh
conda activate snakemake
snakemake -j 6 --use-conda --conda-prefix ~/conda -d test denoise cmp taxonomy ITS
```

A complete denoising run of a dataset on a HPCC may look like this (with SLURM using [this profile](https://github.com/Snakemake-Profiles/slurm#quickstart)):

```sh
outdir=~/path/to/analysis  # must contain a config directory
conda activate snakemake
snakemake -j 10 -c20 \
    --use-conda --conda-prefix ~/conda \
    --profile slurm \
    -d $outdir \
    --rerun-incomplete \
    denoise cmp taxonomy ITS \
    --default-resources runtime=180
```

## Using in R

The R source file [`R/read_input.R`](R/read_input.R) provides code for reading all data from a results directory. See also the small [example analysis](test/R_example/example.md).

## Still not finished...

This repository was created to fit the needs of the author, but lots of things could still be done:

- Enable single-end analyses and other platforms than Illumina
- Integrate more pipelines / clustering methods
- More taxonomy databases
- Better handling of configuration defaults
- Better documentation of configuration options
- Testing deployment on different systems
- ...
