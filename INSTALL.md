# Installation

The Conda package manager is the most important prerequisite. See here for [installation instructions in the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

## Individual pipelines

### QIIME 2

QIIME 2 is installed using Conda according to [these instructions](https://docs.qiime2.org/2022.11/install/native/#install-qiime-2-within-a-conda-environment). Make sure that the correct QIIME version is set in [`config/config.yaml`](config/config.yaml), e.g.:

```yaml
software:
  qiime:
    conda_env: qiime2-2023.7
```

### uvsnake

The software is automatically downloaded from GitHub, just make sure to specify the correct version in the the `software` section of `config.yaml`, e.g.:

```yaml
software:
  usearch:
    binary: usearch  # only needed with program: usearch
  uvsnake:
    snakemake_env: snakemake
    repo:
      tag: v0.1
```

By default, VSEARCH is used, but in case of using USEARCH (`program: usearch` in `config.yaml`). Make sure that it is installed as `usearch` in `$PATH`, or specify the path to USEARCH in the *software* section as shown in the example. USEARCH can be [obtained here](https://www.drive5.com/usearch/download.html)


### Amptk

See [these instructions](https://amptk.readthedocs.io/en/latest/#install) (using Conda). In order to correctly use UNOISE3, make sure to install v1.5.5 or higher.

The UPARSE and UNOISE3 pipelines makes use of the `usearch` section in `config.yaml`, while the DADA2 clustering approach uses the ?`chimera_method` setting in the `dada2` section. In case of using `program: usearch`, the `usearch`/`binary` setting has to be correctly specified. Example:

```yaml
software:
  usearch:
    binary: usearch  # only needed with program: usearch
  amptk:
    conda_env: amptk
```
