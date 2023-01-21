# Installation

The Conda package manager is the most important prerequisite. See here for [installation instructions in the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

Qiime is installed using Conda according to [these instructions](https://docs.qiime2.org/2022.11/install/native/#install-qiime-2-within-a-conda-environment). Make sure that the correct QIIME version is set in [`config/config.yaml`](config/config.yaml).

For Amptk, see [these instructions](https://amptk.readthedocs.io/en/latest/#install) (again using Conda). In order to correctly use UNOISE3, make sure to install v1.5.5 or higher.

All further software is automatically installed when running the pipeline (make sure to specify `--use-conda`).

An exception is *Usearch*. It can be [obtained here](https://www.drive5.com/usearch/download.html) and must be accessible as `usearch` in `$PATH`.
