# Output directories / files

## Example file structure

Usually the results are found in `<my analysis>/results/<workflow name>/data` within analysis directory. `data` is actually a symbolic link to another nested directory and is only present if each workflow in the `results` directory has only a single primer combination (marker) and paired/single-end data is not mixed. Here is a full example result of the `unique_samples`, `denoise`, `taxonomy`, `cmp` and `ITS` rules. The latter runs ITSx, but only for ITS amplicons.

```
📦my_analysis
 ├ 📂 config
 │ ├ 🗋 config.yaml
 │ └ 🗋 taxonomy.yaml
 ├ 📂 input
 │ │ └ ... (contains links or copies of raw sequencing files, used internally)
 ├ 📂 logs
 │ └ ... (log files for every step, helpful for debugging)
 ├ 📂 processing
 │ └ ... (temporary data, sometimes large, remove with 'clean' rule)
 ├ 📂 refdb
 │ └ ... (downloaded reference databases grouped by marker)
 ├ 📂 results
 │ ├ 📂 unoise
 │ │ ├ 🗋 config.yaml
 │ │ ├ 🗋 samples.yaml
 │ │ ├ 🗋 sample_report.tsv
 │ │ 📂 validation
 │ │ ├ ... (MultiQC report(s))
 │ │ ├ 📂 data  [symlink to -> 📂 ITS__ITS3-KYO2...ITS4ngsUni directory]
 │ │ ├ 📂 workflow_usearch_unoise3
 │ │ │ ├ 📂 run1_paired
 │ │ │ │ └ 📂 ITS__ITS3-KYO2...ITS4ngsUni
 │ │ │ │ │ ├ 🗋 denoised.fasta
 │ │ │ │ │ ├ 🗋 denoised.biom
 │ │ │ │ │ ├ 🗋 denoised.hdf5.biom
 │ │ │ │ │ ├ 🗋 denoised_otutab.txt.gz
 │ │ │ │ │ └ 🗋 denoised_search.txt.gz
 │ │ │ │ │ ├ 📂 taxonomy
 │ │ │ │ │ │ ├ 🗋 unite-sintax_usearch-sintax_70.txt.gz
 │ │ │ │ │ │ ├ 🗋 unite-sintax_usearch-sintax_70.biom.gz
 │ │ │ │ │ │ ├ 🗋 unite-sintax_usearch-sintax_70.hdf5.biom.gz
 │ │ │ │ │ │ └ 🗋 unite-qiime_sklearn-sklearn_70.txt.gz
 │ │ │ │ │ │ ├ 🗋 unite-qiime_sklearn-sklearn_70.biom.gz
 │ │ │ │ │ │ ├ 🗋 unite-qiime_sklearn-sklearn_70.hdf5.biom.gz
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
 ├ 📂 unique_samples
 │ └ 📂 demux
 │ └ 🗋 samples.tsv
 │ │ ├ 📂 read_files  (with _1/_2 suffixes to make them unique)
 │ │ │ ├ 🗋 sample1_1_R1.fastq.gz
 │ │ │ ├ 🗋 sample1_2_R2.fastq.gz
 │ │ │ ├ 🗋 sample2_1_R1.fastq.gz
 │ │ │ ├ 🗋 sample2_2_R2.fastq.gz
 │ │ │ └ ... (symlinks to all input files)
```

### Multi-workflow/marker results

 In the case of multiple primer / marker / sequencing layout combinations, the individual results are available using the full path: `results/<workflow name>/workflow_<name>/<run>_<layout>/<marker>__<fwd-primer>...<rev-primer>/` whereby *layout* refers to the paired-end (`paired`) or single-end (not yet implemented) sequencing layout. The following (simplified) directory structure results from running the two workflows named `unoise` and `qiime_dada2`, each with three primer combinations for two markers:

```
📦my_analysis
 ├ 📂 results
 │ ├ 📂 unoise
 │ │ ├ 📂 workflow_usearch_unoise3
 │ │ │ ├ 📂 run1_paired
 │ │ │ │ └ 📂 ITS__ITS3-KYO2...ITS4ngsUni
 │ │ │ │ │  └ ...
 │ │ │ │ └ 📂 ITS__gITS7ngs...ITS4ngsUni
 │ │ │ │ │  └ ...
 │ │ │ │ └ 📂 COI__BF2...BR2
 │ │ │ │ │  └ ...
 │ ├ 📂 qiime_dada2
 │ │ ├ 📂 workflow_qiime_dada2
 │ │ │ ├ 📂 run1_paired
 │ │ │ │ └ ...
 ...
```

A simplified example `config/config.yaml` corresponding to the above setup:

```yaml
input:
  (...)

workflows:
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

- **input**: Contains all input sequence files. Gzip-compressed FASTQ files are directly imported by symlinking (no copying involved). Run pools (`pool_raw: true`) may be generated. This will also be the place for demultiplexed data, once demultiplexing is implemented.
- **unique_samples**: (generated by *unique_samples* rule) Lists all input sample files in one directory and provides `samples.tsv` with useful metadata ([example](../test/unique_samples/demux/samples.tsv)). Samples with identical names (different run/layout combinations) are de-duplicated by adding a numbered suffix. This table is especially useful when submitting sample metadata to public read archives.
- **processing**: Contains all temporary data used by the denoising workflows. In some cases, this can result in gigabytes of data, so be careful to monitor its size. Once everything is finished (deposited in the *results* dir), it can be safely removed with the *clean* command. Downstream target rules such as *taxonomy*, *cmp* and *ITS* don't need this directory to be present.
- **results**: Results directory
  - **sample lists / config files**: Useful to check the input configuration and to obtain sequence metadata (generated by *config*, *collect_input* and *denoise* commands)
    - `results/<workflow>/samples.yaml` lists the sample names and the corresponding input files ([example](../test/results/unoise/samples.yaml))
    - `results/<workflow>/config.yaml` contains the full configuration of each workflow used internally, including configuration defaults not actually visible in `config.yaml` ([example](../test/results/unoise/config.yaml)). This is especially useful to check, whether additional settings specified in the `workflows` definition in `config/config.yaml` correctly override the default settings.
    - `results/<workflow>/sample_report.tsv` lists read statistics for each processing step ([example](../test/results/unoise/sample_report.tsv)).
  - **validation**: Contains QC of every input file and a MultiQC summary (`results/<workflow>/validation/multiqc/multiqc_report.html`). Generated by *quality* and *denoise* command. `multiqc_usearch` further incorporates Cutadapt statistics.
  - **denoising results directories** (*denoise* command): `results/<workflow>/data` in simple cases or `results/<workflow>/workflow_.../<run>_<layout>/<marker>__<fwd-primer>...<rev-primer>` in the case of multiple read layouts and primer combinations. Contents:
    - `denoised.fasta`: Denoised / clustered sequences
    - `denoised_otutab.txt.gz` / `denoised.biom`: OTU table in flat text format (gzip-compressed) or in BIOM format
    - **taxonomy** (*taxonomy* command): The directory containing the taxonomic assignments. The file names are constructed as follows `<ref. database>-<assignment method>-<assignment name>.txt.gz`. In addition, a GZIP-compressed BIOM file is created.
      - *fasta*: subdirectory containing FASTA files annotated with taxonomy
    - **cmp** (*cmp* command): Results of sequence comparisons, which are tab-delimited mapping files with three columns: *query* (ASV/OTU), *target* (from custom sequence database), *percent identity* (edit distance excluding terminal gaps in the alignment, see also `-iddef` USEARCH/VSEARCH option). FASTA files of both the non-matched query and database sequences are also created, as well as a *BAM file*, which allows viewing/extracting the alignments of query and target sequences.
- **logs**: Contains all logfiles by the different commands. Usually it is not necessary to consult them unless there is an error. *Logs* is left in place by the *clean* command.
- **refdb**: Contains taxonomy reference databases. [See here](../workflow/rules/taxonomy.smk) for more information on the structure.
  