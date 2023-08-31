# Test data

The directory `test/gz` contains example files from an Illumina MiSeq run of two uneven (staggered) mock communities (Schlegel et al. 2018). A subset of 5000 reads is stored in `test/gz/run1`. Another 5000 reads from the first uneven community (mock1) are in separate files in the `run2` directory to test run pooling (even though they are not technically from different runs). Due to the reduced sequencing depth, only the more abundant species in the community are present. Still, this dataset is useful for validating the pipelines and their parameter combinations.

## How to analyze

These commands run all the clustering workflows (on a local computer) and compare the results.

```sh
conda activate snakemake

# Run denoising, ITSx, seqence comparisons and taxonomy assignment;
# to make sure that the order of ASVs does not change between runs,
# we use only one core (-c1).
./snakecharmer test denoise ITS cmp taxonomy

# (optional) remove working directories
./snakecharmer test clean

# run a general comparison script (useful for any pipeline comparison)
# (creates Excel file in test/cmp)
scripts/compare_results.sh test
# compare obtained ASVs with mixed gDNA concentrations (see test/mock_cmp/...)
(cd test && scripts/compare_mock.R)
## render the example Rmd (requires pandoc in PATH or RSTUDIO_PANDOC set, here for Ubuntu)
# If this doesn't work, you can still directly run the document in RStudio
export RSTUDIO_PANDOC=/usr/lib/rstudio/resources/app/bin/quarto/bin/tools
Rscript -e "rmarkdown::render('test/R_example/example.Rmd', 'github_document')"

# the following command removes everything (INCLUDING the results/ directory)
./snakecharmer test clean_all
```

The results of the comparison are found in `test/mock_cmp`. The mock community is fairly well represented by the UNOISE, Amptk/UNOISE and Amptk/Dada2 pipelines, while QIIME2 looses many reads due to quality filtering. Using this pipeline with variable-length ITS reads seems difficult.

![mock comparison](mock_cmp/ITS__ITS3-KYO2...ITS4/mock.png)


## Validation using example workflows

In order to further carefully validate this software, the test data was further analyzed using example scripts from the online documentation of the different tools, currently:

* [USEARCH pipeline](https://www.drive5.com/usearch/manual/ex_miseq_its.html) for MiSeq 2x300 fungal ITS
* [VSEARCH "alternative" pipeline](https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline/c4859786f05bba35d8c306de4a3d64fea40d9dbf) slightly modified to use UNOISE3 following the [this description](https://github.com/torognes/vsearch/pull/283). The VSEARCH "alternative" pipeline contains an extra step of read mapping against the OTUs to obtain the count table, using quality filtered reads in this case. The workflow from this repository maps the raw/unfiltered reads instead (with a 97% identity threshold), [as recommended by the USEARCH author](https://www.drive5.com/usearch/manual/cmd_otutab.html). In the future, this should be configurable.

The following script runs the "simple" workflows and compares the results with the outcomes of our test pipeline:

```sh
# if running for the first time, do this first:
# conda env create -f scripts/simple/uvsearch_env.yaml
scripts/simple/compare.sh
```

## Reference

Schlegel, M., Queloz, V., and Sieber, T. N. (2018). The endophytic mycobiome of European ash and sycamore maple leaves – geographic patterns, host specificity and influence of ash dieback. *Frontiers in Microbiology* 9. doi: 10.3389/fmicb.2018.02345.
