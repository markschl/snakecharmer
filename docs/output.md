# Output directories / files

## Example file structure

Usually the results are found in `<my analysis>/results/<workflow name>/data` within analysis directory. `data` is actually a symbolic link to another nested directory and is only present if each workflow in the `results` directory has only a single primer combination (marker) and paired/single-end data is not mixed. Here is a full example result of the `denoise`, `taxonomy`, `cmp` and `ITS` rules. The latter runs ITSx, but only for ITS amplicons.

```
📦my_analysis
 ├ 📂 config
 │ ├ 🗋 config.yaml
 │ └ 🗋 taxonomy.yaml
 ├ 📂 input
 │ ├ 📂 grouped
 │ │ └ ... (contains links or copies of raw sequencing files, used internally)
 │ └ 📂 unique_samples
 │ │ ├ 🗋 sample1_R1.fastq.gz
 │ │ ├ 🗋 sample1_R2.fastq.gz
 │ │ ├ 🗋 sample2_R1.fastq.gz
 │ │ ├ 🗋 sample2_R2.fastq.gz
 │ │ └ ... (symlinks to all input files)
 ├ 📂 logs
 │ └ ... (log files for every step, helpful for debugging)
 ├ 📂 processing
 │ └ ... (temporary data, sometimes large, remove with 'clean' rule)
 ├ 📂 refdb
 │ └ ... (downloaded reference databases grouped by marker)
 ├ 📂 results
 │ ├ 📂 unoise
 │ │ ├ 📂 data  [symlink to -> 📂 paired directory]
 │ │ ├ 📂 pipeline_usearch_unoise3
 │ │ │ ├ 📂 ITS__ITS3-KYO2...ITS4ngsUni
 │ │ │ │ └ 📂 paired
 │ │ │ │ │ ├ 🗋 denoised.fasta
 │ │ │ │ │ ├ 🗋 denoised.biom
 │ │ │ │ │ ├ 🗋 denoised_otutab.txt.gz
 │ │ │ │ │ └ 🗋 denoised_search.txt.gz
 │ │ │ │ │ ├ 📂 taxonomy
 │ │ │ │ │ │ ├ 🗋 unite-sintax_usearch-sintax_70.txt.gz
 │ │ │ │ │ │ ├ 🗋 unite-sintax_usearch-sintax_70.biom.gz
 │ │ │ │ │ │ └ 🗋 unite-qiime_sklearn-sklearn_70.txt.gz
 │ │ │ │ │ │ ├ 🗋 unite-qiime_sklearn-sklearn_70.biom.gz
 │ │ │ │ │ │ ├ 📂 fasta
 │ │ │ │ │ │ │ ├ 🗋 unite-sintax_usearch-sintax_70.fasta.gz
 │ │ │ │ │ │ │ └ 🗋 unite-qiime_sklearn-sklearn_70.fasta.gz
 │ │ │ │ │ │ ├ 📂 sintax
 │ │ │ │ │ │ │ └ ... (original SINTAX output)
 │ │ │ │ │ ├ 📂 cmp
 │ │ │ │ │ │ ├ 🗋 my_file_comparison.txt
 │ │ │ │ │ │ ├ 🗋 my_file_comparison.bam
 │ │ │ │ │ │ ├ 🗋 my_file_comparison.bam.bai
 │ │ │ │ │ │ ├ 🗋 my_file_comparison_denoised_notmatched.fasta.gz
 │ │ │ │ │ │ └ 🗋 my_file_comparison_notmatched.fasta.gz
 │ │ │ │ │ ├ 📂 [ITSx]
 │ │ │ │ │ │ ├ 🗋 out.positions.txt
 │ │ │ │ │ │ └ ... (ITSx output)
 │ │ │ ├ 📂 _validation
 │ │ │ │ ├ 📂 multiqc
 │ │ │ │ │ └ 🗋 multiqc_report.html  [including primer-trimming report]
 │ │ │ │ └ 🗋 sample_report.tsv
 │ ├ 📂 _validation
 │ │ └ ... (FastQC / MultiQC results)
 │ ├ 🗋 samples.tsv
 │ └ 🗋 samples.yaml
```

### Multi-workflow/marker results

 In the case of multiple primer / marker / sequencing strategy combinations, the individual results are available using the full path: `results/<workflow name>/pipeline_<name>/<marker>__<fwd-primer>...<rev-primer>/<strategy>/` whereby *strategy* refers to the paired-end (`paired`) or single-end (not yet implemented) sequencing strategy. The following (simplified) directory structure results from running the two workflows named `unoise` and `qiime_dada2`, each with three primer combinations for two markers:

```
📦my_analysis
 ├ 📂 results
 │ ├ 📂 _validation
 │ │ └ ... (FastQC / MultiQC results)
 │ ├ 📂 unoise
 │ │ ├ 📂 pipeline_usearch_unoise3
 │ │ │ ├ 📂 ITS__ITS3-KYO2...ITS4ngsUni
 │ │ │ │ └ 📂 paired
 │ │ │ │ │  └ ...
 │ │ │ ├ 📂 ITS__gITS7ngs...ITS4ngsUni
 │ │ │ │ └ 📂 paired
 │ │ │ │ │  └ ...
 │ │ │ ├ 📂 COI__BF2...BR2
 │ │ │ │ └ 📂 paired
 │ │ │ │ │  └ ...
 │ │ └ 🗋 config.yaml
 │ ├ 📂 qiime_dada2
 │ │ ├ 📂 pipeline_qiime_dada2
 │ │ │ ├ 📂 ITS__ITS3-KYO2...ITS4ngsUni (...)
 │ │ │ ├ 📂 ITS__gITS7ngs...ITS4ngsUni (...)
 │ │ │ ├ 📂 COI__BF2...BR2 (...)
 │ │ └ 🗋 config.yaml
 │ ├ 🗋 samples.tsv
 │ └ 🗋 samples.yaml
```

A simplified example `config/config.yaml` corresponding to the above setup:

```yaml
input:
  (...)

pipelines:
  unoise:
    cluster: usearch_unoise3
    taxonomy: default
  qiime_dada2:
    cluster: qiime_dada2
    taxonomy: default

compare:
  my_file_comparison:
    file: my_sequences.fasta
    ident_threshold: 0.9

primers:
  ITS:
    forward: 
      - ITS3-KYO2: GGGATGAAGAACGYAGYRAA
      - gITS7ngs: GTGARTCATCRARTYTTTG
    reverse:
      - ITS4ngsUni: CCTSCSCTTANTDATATGC
    combinations: default  # both forward primers will paired with the one reverse primer
  COI:
    forward:
      -BF2: GCHCCHGAYATRGCHTTYCC
    reverse:
      -BR2: TCDGGRTGNCCRAARAAYCA
    combinations: default
  (...)

taxonomy_dbs:
  ITS:
    unite:
      db: unite_eukarya_all
      defined: species
  COI:
    midori:
      db: Midori_COI
      defined: species

taxonomy_methods:
  sintax_70:
    method: sintax_usearch
    confidence: 0.7
  sklearn_70:
    method: qiime_sklearn
    confidence: 0.7

(...)
```

## Detailed description of directories/files

- **input**: Contains all input sequence files (FASTQ, generated by *collect_input* and many other commands)
  - **input/grouped**: hierarchical grouping of the files, which is then used as input for the pipelines.
  - **input/unique_samples**: lists all input sample files in one directory. Samples with identical names (but from different directories) are de-duplicated by adding a numbered suffix.
- **processing**: Contains all temporary data used by the denoising pipelines. In some cases, this can result in gigabytes of data, so be careful to monitor its size. Once everything is finished (deposited in the *results* dir), it can be safely removed with the *clean* command. Downstream target rules such as *txonomy*, *cmp* and *ITS* don't need this directory to be present.
- **results**: Results directory
  - **sample lists / config files**: Useful to check the input configuration and to obtain sequence metadata (generated by *config*, *collect_input* and *denoise* commands)
    - `samples.yaml` lists the sample names and the corresponding input files ([example](../test/results/samples.yaml))
    - The tab-separated `results/samples.tsv` file also lists all unique samples and their corresponding read files in tabular format ([example](../test/results/samples.tsv)). This table is especially useful when submitting sample metadata to public SRA archives. For files from different directories with the same sample name, a numbered suffix is added to the sample to make the names unique.
    - The files `<pipeline>/config.yaml` contain the final configuration of each pipeline used internally (slightly modified, structure may change with later versions) ([example](../test/results/unoise/config.yaml)). This is especially useful to check, whether additional settings specified in the `pipelines` definition in `config/config.yaml` correctly override the default settings.
  - **results/_validation**: Contains QC of every input file and a MultiQC summary (`results/_validation/multiqc/multiqc_report.html`). Generated by *quality* and *denoise* command
  - **denoising results directories** (*denoise* command): `results/<pipeline>/data` in simple cases or `results/<pipeline>/pipeline_.../<marker>__<fwd-primer>...<rev-primer>/<strategy>` in the case of multiple primer combinations and sequencing strategies. Contents:
    - `denoised.fasta`: Denoised / clustered sequences
    - `denoised_otutab.txt.gz` / `denoised.biom`: OTU table in flat text format (gzip-compressed) or in BIOM format
    - **taxonomy** (*taxonomy* command): The directory containing the taxonomic assignments. The file names are constructed as follows `<ref. database>-<assignment method>-<assignment name>.txt.gz`. In addition, a GZIP-compressed BIOM file is created.
      - *fasta*: subdirectory containing FASTA files annotated with taxonomy
    - **cmp** (*cmp* command): Results of sequence comparisons, which are tab-delimited mapping files with three columns: *query* (ASV/OTU), *target* (from custom sequence database), *percent identity* (edit distance excluding terminal gaps in the alignment, see also `-iddef` USEARCH/VSEARCH option). FASTA files of both the non-matched query and database sequences are also created, as well as a *BAM file*, which allows viewing/extracting the alignments of query and target sequences.
  - **sample report** (*denoise* command): A sample report with read numbers retained in each step is placed in `results/<pipeline>/pipeline_.../_validation`, along with another MultiQC report including Cutadapt statistics (currently only UNOISE pipeline).
- **logs**: Contains all logfiles by the different commands. Usually it is not necessary to consult them unless there is an error. *Logs* is left in place by the *clean* command.
- **refdb**: Contains Zstandard-compressed taxonomy reference databases, including trained reference sets.