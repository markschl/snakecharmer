# Target rules ("commands") and output files

## List of rules

### Initial checks

- **config**: Allows checking the workflow configuration, creates the files `results/samples.yaml` and `<workflow>/config.yaml`
- **samples**: Parses the *input* configuration in `config/config.yaml` and creates links (symlinks) of the input files to the `input` directory. Lists all samples in `results/<workflow>/samples.yaml`.
- **unique_samples**: Creates a directory `unique_samples/demux/read_files` with symlinks to all input files with unique names (numbered suffix appended in case of multiple sequencing runs with same samples). The file `unique_samples/demux/read_files` contains the necessary metadata (including MD5 hashes) e.g. for upload to public read archives.
- **quality**: Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) and [MultiQC](https://multiqc.info) (output in `results/<workflow>/validation`). This can be done in a first step before deciding on the quality filtering/trimming options in `config/config.yaml`.

### Denoising (clustering)

- **denoise**: Does the denoising (clustering) for all workflows, as specified in `config/config.yaml`. The results are normally found in `results/<workflow>/data`, unless the input files contain both paired and single-end (different "sequencing strategies"), or multiple primer combinations are specified. In this case, the output path is `results/<workflow>/workflow_.../<run>_<layout>/<primers>`.

### Analysis of denoised sequences

- **taxonomy**: Assigns taxonomic labels to denoised/clustered sequences. Several methods can be specified in `config/config.yaml`. The output files are in the `taxonomy` subdirectory of the denoising results directory.
- **cmp**: Runs sequence comparisons with files listed in `config/config.yaml` (`compare` key). The comparisons are done with VSEARCH, using the denoised FASTA as query and the given sequence files as database. The result is a tab-delimited mapping file stored in the `cmp` subdirectory of the denoising results directory (along with a few other files).
- **ITS**: Runs [ITSx](https://microbiology.se/software/itsx) to recognize rDNA domains (or parts of them) and locate the internal transcribed spacer (ITS) regions. The information from the positions file (`ITSx/out.positions.txt`) is used to distinguish "true" ITS sequences from possible unspecific amplification (see also [example analysis](test/R_example/example.md#read-data)).

### Cleanup

- **clean**: Removes the working directories `input`, `processing` and `logs`. The `results` and `unique_samples` directories are retained.
- **clean_tax**, **clean_cmp**, **clean_itsx**: Removes the `taxonomy`, `cmp` or `ITSx` directories in all results directories.
- **clean_all**: Cleans up everything (including the workflow output). This should only be used to **completely remove** all output from a target directory and obtain a clean workspace.
