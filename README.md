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

As visible in the rulegraph above, there are a few Snakemake target rules, which can be run independently (even though they also partly depend on each other):

- `config`: Allows checking the configuration by creating the files `results/samples.yaml` and `<pipeline>/config.yaml`. The samples file (`samples.yaml`) lists the sample names and the corresponding input files. This is useful to check if the `input` settings in `config/config.yaml` are specified correctly. The `config.yaml` files contain the final configuration of each pipeline used internally (slightly modified, structure may change with later versions). This is especially useful to check, whether additional settings specified in the `pipelines` definition in `config/config.yaml` correctly override the default settings.
- `quality`: Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) and [MultiQC](https://multiqc.info) (output in `results/_validation`). This can be done in a first step before deciding on the quality filtering/trimming options in `config/config.yaml`.
- `denoise`: Does the denoising (clustering) for all pipelines, as specified in `config/config.yaml`. The results are normally found in `results/<pipeline>/data`, unless the input files contain both paired and single-end (different "sequencing strategies"), or multiple primer combinations are specified. In this case, the output path is `results/<pipeline>/pipeline_.../<primers>/<strategy>`. The denoised/clustered sequences are in `denoised.fasta`, the OTU table is stored in flat text format (`denoised_otutab.txt.gz`) and BIOM (`denoised.biom`). In addition, a sample report with read numbers retained in each step is placed in `results/<pipeline>/pipeline_.../_validation`, along with another MultiQC report including Cutadapt statistics (currently only UNOISE pipeline).
- `taxonomy`: Assigns taxonomic labels to denoised/clustered sequences. Several methods can be specified in `config/config.yaml`. The output files are in the `taxonomy` subdirectory of the directory containing the FASTA and OTU tables. The file names are as follows `<ref. database>-<assignment method>-<assignment name>.txt.gz`. In addition, a GZIP-compressed BIOM file is created.
- `cmp`: Runs sequence comparisons with files listed in `config/config.yaml` (`compare` key). The comparisons are done with VSEARCH, using the denoised FASTA as query and the given sequence files as database. The result is a tab-delimited mapping file (with corresponding name) with three columns: query (ASV/OTU), target (from custom sequence database), percent identity (edit distance excluding terminal gaps in the alignment, see also --iddef option). FASTA files of both the non-matched query and database sequences are also created, as well as a BAM file, which allows viewing/extracting the alignments of query and target sequences.
- `ITS`: Runs [ITSx](https://microbiology.se/software/itsx) to recognize rDNA domains (or parts of them) and locate the internal transcribed spacer (ITS) regions. The information from the positions file (`ITSx/out.positions.txt`) is used to distinguish "true" ITS sequences from possible unspecific amplification (see also `R/read_amplicon.R` file).

Finally, there are commands for cleaning up:

- `clean`: Removes the working directories `input` and `processing`. The `results` and `logs` directories are retained.
- `clean_all`: Cleans up everything (including the pipeline output), except for taxonomic reference databases. This is only used to completely remove all output from a target directory.
- `clean_tax`: Removes the `refdb` directory

### On a local computer

The following command runs the test pipeline (FASTQ files from sequencing fungal mock comunities in the [`test` directory](test/) directory) on a local computer.

```sh
conda activate snakemake
snakemake -j 6 --use-conda --conda-prefix ~/conda -d test denoise cmp taxonomy ITS
```

Note that the `~/conda` directory is used for the installation of all additional necessary software. This allows reusing the installed software across analyzes of different datasets.

### On computing cluster

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

The R source file [`R/read_amplicon.R`](R/read_amplicon.R) provides code for reading all data from a results directory. See also the small [example analysis](test/R_example/example.md).


## Comparison of denoising/clustering pipelines

There is a separate bash script `scripts/compare_results.sh`, which creates an Excel file comparing the number of reads assigned to each 98% cluster by each pipeline. A separate workbook is created for each sample. The script requires VSEARCH, as well as R with the following packages: `ggplot2`, `tidyverse`, `data.table` and `openxlsx`.


## Still not finished...

This repository was created to fit the needs of the author, but lots of things could still be done:

- Enable single-end analyses and other platforms than Illumina
- Integrate more pipelines / clustering methods
- More taxonomy databases
- Better handling of configuration defaults
- Better documentation of configuration options
- Testing deployment on different systems
- ...
