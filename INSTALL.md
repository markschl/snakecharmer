# Installation

The Conda package manager is the most important prerequisite. See here for [installation instructions in the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

## Individual pipelines

* **QIIME 2** is installed using Conda according to [these instructions](https://docs.qiime2.org/2022.11/install/native/#install-qiime-2-within-a-conda-environment). Make sure that the correct QIIME version is set in [`config/config.yaml`](config/config.yaml).
* **USEARCH/VSEARCH-based pipeline**: By default, VSEARCH is used, but in case of using USEARCH (`program: usearch`), make sure that it is installed as `usearch` in `$PATH`. Alternatively, specify the path to USEARCH in the *software* section of `config.yaml`. USEACH can be [obtained here](https://www.drive5.com/usearch/download.html).
* **Amptk**: see [these instructions](https://amptk.readthedocs.io/en/latest/#install) (again using Conda). In order to correctly use UNOISE3, make sure to install v1.5.5 or higher.
