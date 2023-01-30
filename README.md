# Amplicon data processing workflows for diversity analysis

This software makes use of the workflow management system [Snakemake](https://snakemake.github.io/) to build and integrate different amplicon workflows and allows for the flexible processing of data from different markers.

**Features**

- Denoising/clustering of raw amplicon sequencing data, taxonomic assignments and further sequence comparisons and marker-specific processing (such as ITS extraction)
- Comparison of different pipelines or variations of the same pipeline using different sets of parameters based on a flexible [`configuration system`](config/config.yaml). Downstream analyses and comparisons are easily possible due to a common file structure of the output directories.
- Simultaneous processing of multi-marker amplicons generated using different primer sets
- Multiple taxonomic assignment methods can be applied to each denoised dataset using marker-specific reference databases; currently implemented: [UNITE](https://unite.ut.ee) for Eukaryote ITS and [Midori](http://www.reference-midori.info) for mitochondrial markers (more may follow)

*Note:* To this date, a few workflows for paired-end Illumina data and taxonomic assignments for Eukaryote lineages have been implemented, but the software will be extended further (see [below](#further-steps)).

No snakes 🐍 were harmed in the process of creating this software

**Non-features**

The software is mainly intended for the basic tasks of obtaining an OTU table, taxonomic assignments, sequence searches and basic validation of different workflows, even though the integrated pipelines often offer additional comprehensive downstream analyses. However, the output files are in commonly used formats (such as BIOM) and can serve as input for many analysis toolkits. There is also an [import script for R](#analyzing-in-r).

**Integrated pipelines**

- [USEARCH](https://www.drive5.com/usearch/manual)/[VSEARCH](https://github.com/torognes/vsearch)-based amplicon pipeline using UNOISE3 for obtaining ASVs
- [QIIME2](https://qiime2.org) (currently with DADA2 denoising)
- [Amptk](https://github.com/nextgenusfs/amptk) (UNOISE3 and DADA2)
- ... (more to follow)

**Validation**

Validation is done using amplicon data from a fungal mock community ([details in `test` directory](test/README.md)) and a basic comparison of the different pipelines can be done [using a script](#comparison-of-denoisingclustering-pipelines).

## Installing

The pipeline makes use of the [Conda package manager](https://conda.io), the installation is thus pretty straightforward ([see instructions here](INSTALL.md)).

## Configuring

The easiest is to copy the contents of the [config](config/) directory (or [test/config](test/config)) into the future analysis directory and then modify the files according to your needs. [`config.yaml`](config/config.yaml) contains all settings, while the taxonomic databases can be configured in [`taxonomy.yaml`](config/taxonomy.yaml). For the time being, please refer to the comments in these config files.

## Running

There are a few Snakemake target rules, which can be run independently (even though they also partly depend on each other). A complete list of commands is [found here](docs/rules.md).

### On a local computer

The following command runs the test pipeline (FASTQ files from sequencing fungal mock comunities in the [`test` directory](test/), specified with `-d test`) using 6 cores on a local computer. The  [target rules](docs/rules.md) to be run are `denoise`, `cmp`, `taxonomy` and `ITS`.

```sh
conda activate snakemake
snakemake -c6 --use-conda --conda-prefix ~/conda -d test denoise cmp taxonomy ITS
```

Note that the `~/conda` directory is used for the installation of all additional necessary software (`--conda-prefix` argument). This allows reusing the installed software across analyzes of different datasets.

### Output

After running, a few additional directories will have appeared next to `config`. The most important one is the `results` directory, which roughly has the following structure ([more details here](docs/output.md)):

```
📦<my_analysis>
 ├─ 📂 config
 │  ├─ 🗋 config.yaml
 │  └─ 🗋 taxonomy.yaml
 │  (...)
 ├─ 📂 results
 │  ├─ 📂 <workflow name>
 │  │  ├─ 📂 data
 │  │  │  ├─ 🗋 denoised.fasta
 │  │  │  ├─ 🗋 denoised_otutab.txt.gz
 │  │  │  ├─ 🗋 denoised.biom
 │  │  │  ├─ 🗋 denoised_search.txt.gz
 │  │  │  ├─ 📂 taxonomy
 │  │  │  │  ├─ 🗋 <database>-<method>-<name>..txt.gz
 │  │  │  │  ├─ 🗋 <database>-<method>-<name>.biom.gz
 │  │  │  │  │  (...)
 │  │  │  ├─ 📂 cmp
 │  │  │  │  ├─ 🗋 <my_seq_comparison>.txt
 │  │  │  │  ├─ 🗋 <my_seq_comparison>_notmatched.fasta.gz
 │  │  │  │  ├─ 🗋 <my_seq_comparison>_denoised_notmatched.fasta.gz
 │  │  │  │  │  (...)
 │  │  │  ├─ 📂 [ITSx]
 │  │  │  │  ├─ 🗋 out.positions.txt
 │  │  │  │  └─ (...)
 │  │  ├─ 📂 pipeline_<type>
 │  │  │  ├─ 📂 <marker>__<fwd-primer>...<rev-primer>
 │  │  │  │  └─ 📂 <single/paired>
 │  │  │  │     └─ (... same as in *data* directory, only relevant with 
 │  │  │  │             multi-marker/workflow setups)
```

Whether the output are ASVs/ESVs or OTUs from a fixed threshold clustering (not yet implemented), the resulting FASTA file is always called `denoised.fasta`. The sample/OTU count matrix is returned both in the traditional tabular format (`denoised_otutab.txt.gz`) and [BIOM v1](https://biom-format.org/documentation/biom_format.html). The taxonomic annotations are named by taxonomy database and assignment method (multiple combinations possible) and returned in a QIIME-style tabular format as well as the BIOM format. Furthermore, there can be results of sequence comparisons (`cmp`) or marker-specific data such as ITSx results.

With the most simple scenario (one workflow with one primer combination), the relevant results directory is `<my_analysis>/results/<pipeline_name>/data`. With multi-workflow/marker setups, the `data` directory will not be present, and the individual workflow results are placed in the nested directories.

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

For the USEARCH-based pipeline, it is recommended to use `--group-components` to limit the number of submitted jobs. Let's assume that 300 paired-end samples should be analyzed. Specifying `--group-components sample=50` will make sure that paired-end read merging, primer trimming and quality filtering is done in batches of 50 samples, and thus only six sample processing jobs have to be submitted to the cluster instead of 300. The same can be done for QC jobs on individual read files by specifying `qc=100` (again to obtain six jobs). After pre-processing, the sequences are combined for denoising, which can only run on a single node. However, if many workflows should be evaluated simultaneously, component grouping can also be used (e.g. `--group-components denoise=4`).

Other pipelines (currently QIIME and Amptk) don't support processing samples in parallel on different nodes. For these, the pre-processing steps belong to another group called `prepare`. Furthermore, the following groups exist: `otutab` (OTU table construction for USEARCH pipelines), `taxonomy` (taxonomy assignment), `ITS` (ITSx) and `cmp` (sequence comparisons). It is also possible to [assign single rules to custom groups](https://snakemake.readthedocs.io/en/v7.19.1/executing/grouping.html#job-grouping).

## Analyzing in R

The R source file [`R/read_amplicon.R`](R/read_amplicon.R) provides code for reading all data from a results directory. See also the small [example analysis](test/R_example/example.md).


## Comparison of denoising/clustering pipelines

There is a separate bash script `scripts/compare_results.sh`, which creates an Excel file comparing the number of reads assigned to 98% clusters of the already denoised sequences by each pipeline. A separate workbook is created for each sample. The script requires VSEARCH, as well as R with the following packages: `ggplot2`, `tidyverse`, `data.table` and `openxlsx`.


## Further steps...

A list of possible next steps includes:

- Integrate more pipelines / clustering methods and taxonomy databases
- Integrate other platforms than Illumina and allow simultaneous analysis of multi-platform data
- Allow for separate analysis of sequencing runs and result merging
- Offer more ways of comparing and validating pipelines and generally improve user experience
- Testing deployment on different systems
- Improve configuration of job resources (memory, CPUs)
- The USEARCH pipeline may be moved into an extra repository to be used independently
- ...

