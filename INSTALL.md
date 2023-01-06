# Installation

The Conda package manager is the most important prerequisite. See here for [installation instructions in the Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

Qiime is installed using Conda according to [these instructions](https://docs.qiime2.org/2022.8/install/native/#install-qiime-2-within-a-conda-environment). Make sure that the correct QIIME version is set in [`config/config.yaml`](config/config.yaml).

For Amptk, see [these instructions](https://amptk.readthedocs.io/en/latest/#install) (again using Conda). In addition, there is currently a [bug](https://github.com/nextgenusfs/amptk/issues/96) with UNOISE3 denoising. The code can be fixed like this:

```sh
conda activate amptk
sed -i "s/--derep_fulllength', filter_out, '--relabel', 'Read_', '--sizeout', '--output/--fastx_uniques', filter_out, '--relabel', 'Read_', '--sizeout', '--fastaout/g" "$CONDA_PREFIX/lib/python3.10/site-packages/amptk/unoise3.py"
```

All further software is automatically installed when running the pipeline (make sure to specify `--use-conda`).

An exception is *Usearch*. It can be [obtained here](https://www.drive5.com/usearch/download.html) and must be accessible as `usearch` in `$PATH`.

