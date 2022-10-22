# Analysis of environmental diversity using amplicon sequencing

This pipeline system uses [Snakemake](https://snakemake.github.io/) to create and integrate different amplicon pipelines. Their results can be easily compared due to a common file structure of the output directories, which may serve as input e.g. for R-based analyses.

*Note:* The repository was created to fit the needs of the author and is tailored for the analysis of paired-end Illumina data from fungi or arthropods, but it may be extended further (see [below](#still-not-finished)).

*Pipelines*

- [USEARCH](https://www.drive5.com/usearch/manual)/[VSEARCH](https://github.com/torognes/vsearch)-based amplicon pipeline using UNOISE3 for obtaining ASVs (implemented here from scratch with Snakemake)
- [QIIME2](https://qiime2.org) (DADA2 denoising)
- [Amptk](https://github.com/nextgenusfs/amptk) (UNOISE3 and DADA2)
- (more may follow)

*Features*

- Running different pipelines or variations of the same pipeline using different sets of parameters, which are defined in a [`configuration file`](config/config.yaml).
- Simultaneous processing of amplicons generated using different primer sets
- Multiple taxonomic assignment methods can be applied to the denoised dataset using marker-specific reference databases; currently implemented: [UNITE](https://unite.ut.ee) for Eukaryote ITS and [Midori](http://www.reference-midori.info) for mitochondrial markers (more may follow)

The pipeline is tested using amplicon data from a fungal mock community ([`test` directory](test/)) and the results of the different pipelines can be compared [using a script](#comparison-of-denoisingclustering-pipelines).

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

The easiest is to copy the contents of the [config](config/) directory into the future analysis directory and then modify the files according to your needs. [`config.yaml`](config/config.yaml) contains all settings, while the taxonomic databases can be configured in [`taxonomy.yaml`](config/taxonomy.yaml). Please refer to the comments in these config files.

## Running

As visible in the rulegraph above, there are a few Snakemake target rules, which can be run independently (even though they also partly depend on each other). A complete list of commands is [found here](Commands.md).

### On a local computer

The following command runs the test pipeline (FASTQ files from sequencing fungal mock comunities in the [`test` directory](test/)) on a local computer.

```sh
conda activate snakemake
snakemake -j 6 --use-conda --conda-prefix ~/conda -d test denoise cmp taxonomy ITS
```

Note that the `~/conda` directory is used for the installation of all additional necessary software. This allows reusing the installed software across analyzes of different datasets.

### On a computer cluster

A complete denoising run of a dataset on a HPCC may look like this (using SLURM with [this profile](https://github.com/Snakemake-Profiles/slurm#quickstart)):

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


## Analyzing in R

The R source file [`R/read_amplicon.R`](R/read_amplicon.R) provides code for reading all data from a results directory. See also the small [example analysis](test/R_example/example.md).


## Comparison of denoising/clustering pipelines

There is a separate bash script `scripts/compare_results.sh`, which creates an Excel file comparing the number of reads assigned to each 98% cluster by each pipeline. A separate workbook is created for each sample. The script requires VSEARCH, as well as R with the following packages: `ggplot2`, `tidyverse`, `data.table` and `openxlsx`.


## Still not finished...

Lots of things could still be done. A list of possible next steps includes:

- Enable single-end analyses and other platforms than Illumina
- Integrate more pipelines / clustering methods
- Integrate more taxonomy databases
- Better handling of configuration defaults
- Better documentation of configuration options, clear error messages in case of problems with the config files.
- Testing deployment on different systems
- ...
