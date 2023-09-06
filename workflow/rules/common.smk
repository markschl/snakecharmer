from os.path import *
import shutil
import lib
import os
from os.path import dirname
import csv
import re
from collections import OrderedDict, defaultdict
from glob import glob
import copy

from snakemake.workflow import srcdir

# Set environment variable of workflow root dir to allow post-deploy
# scripts to run other scripts stored in that directory
os.environ['PIPELINE_DIR'] = dirname(dirname(dirname(srcdir('.'))))


localrules:
    collect_sample_lists,
    make_pooling_list,
    dump_config,
    link_data_dir,
    combine_sample_reports,
    collect_unique_files


ruleorder:
    # If pool_raw is true, an artificial run called 'pool' will be created,
    # resulting in an ambiguous situation, where run pooling needs to have priority
    pool_runs_raw > collect_sample_files


#### Configuration ####

cfg = lib.Config(config)

# commonly used
run_config = "input/sample_config/{technology}/{layout}/{run}"


#### Helper functions related to rules ####


def with_default(_config, group, setting):
    value = _config[group].get(setting)
    out = _config["defaults"][setting] if value is None else value
    assert out is not None
    return out


def expand_input_files(path=None, **wildcards):
    """
    Lists Fastq sample files used for pipeline input, by checking both the
    sample table and the files actually present to be sure that everything is 
    consistent
    """
    wildcards.update(cfg.get_run_data(wildcards["workflow"], wildcards["run"], wildcards["layout"]))
    tab = checkpoints.final_sample_tab.get(**wildcards).output.tab
    # sample_dir = checkpoints.link_input.get(**wildcards).output.sample_dir
    sample_dir = rules.link_input.output.sample_dir.format(**wildcards) + "/nested"
    full_path = sample_dir + "/{sample}_R{read}.fastq.gz"
    w = glob_wildcards(full_path)
    d = cfg.read_samples(tab, **wildcards)
    assert sorted(set(w.sample)) == sorted(d["sample"]) and \
        sorted(set(w.read)) == sorted(d["read"]), (
        "Sample tab does not align with actual files, "
        "try deleting processing/<workflow>/input and re-run")
    return expand(full_path if path is None else path, **wildcards, **d)


def expand_runs(path, workflows=cfg.workflows, pooled=True, **param):
    for workflow in workflows:
        for run_data in cfg.get_runs(workflow, pooled=pooled):
            yield from expand(
                path,
                workflow=workflow,
                cluster=cfg.workflows[workflow]["cluster"],
                **run_data,
                **param
            )


def run_results(sub_path="", workflows=cfg.workflows, **param):
    return expand_runs("results/{workflow}/workflow_{cluster}/{run}_{layout}{sub_path}", workflows=workflows, sub_path=sub_path, **param)


# assists in listing results files
# requires sub-path in the results dir, or a function that accepts all workflow settings and returns a sub-path
# TODO: redundancy with run_results
def result_paths(sub_path="", workflows=cfg.workflows, pooled=True, allow_missing=True):
    for workflow, p in cfg.workflows.items():
        for marker, primer_comb in cfg.primer_combinations.items():
            for r in cfg.get_runs(workflow, pooled=pooled):
                yield from expand(
                    "results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}{sub_path}",
                    workflow=workflow,
                    cluster=p["cluster"],
                    marker=marker,
                    primers=primer_comb,
                    sub_path=sub_path(marker, primer_comb, p) if callable(sub_path) else sub_path,
                    allow_missing=allow_missing,
                    **r
                )


#### Prepare / configure ####

# Deposits the sample file lists ("manifest files") for every run and run layout
# in a standardized directory structure, so they can be used by Snakemake.
# Their content (the samples) will only be known to Snakemake later after the
# link_input checkpoint
# Most importantly, reserved characters in sample names are replaced
rule collect_sample_lists:
    params:
        run_meta = lambda _: list(cfg.get_runs(pooled=False)),
        reserved_chars = ".- ",  # TODO: hardcoded, should it be configurable?
        path_template = lambda _: run_config + "/samples.{ext}",
    output:
        sample_files=sorted(chain(*[
            expand(
                run_config + "/samples.{ext}",
                ext=["tsv", "yaml"],
                **d
            )
            for d in cfg.get_runs(pooled=False)
        ])),
    log:
        "logs/prepare/collect_sample_lists.log",
    script:
        "../scripts/collect_sample_lists.py"


rule dump_config:
    params:
        config=lambda wildcards: cfg.workflows[wildcards.workflow],
    output:
        "results/{workflow}/config.yaml",
    log:
        "logs/prepare/{workflow}/dump_config.log",
    script:
        "../scripts/dump_config.py"


# checks and reverse complements primer sequences and calculates the consensus
# of primer mixtures
rule prepare_primers:
    params:
        primers=cfg.primers,
    output:
        yaml="processing/primers/primers.yaml",
    log:
        "logs/prepare/prepare_primers.log",
    conda:
        "envs/consensus.yaml"
    script:
        "../scripts/prepare_primers.py"


#### Demultiplexing ####

# # TODO: nonfunctional stub

# # adds non-demultiplexed R1 [and R2] files to the input dir as symlinks, given
# # input options
# rule collect_raw:
#     input:
#         indexes=run_config + "/indexes.tsv"
#     output:
#         fq=expand(
#             "input/{{technology}}/{{layout}}/{{run}}/raw/R{read}.fastq.gz",
#             read=[1, 2]
#         ),
#     run:
#         pass


# # example demultiplexing rule
# rule demux_prog1:
#     input:
#         fq=expand(
#             "input/{{technology}}/{{layout}}/{{run}}/raw/R{read}.fastq.gz",
#             read=[1, 2]
#         ),
#         indexes=run_config + "/indexes.tsv"
#     output:
#         fq=directory("input/{technology}/{layout}/demux_prog1_{method}/{run}"),
#     run:
#         pass


#### Sample handling ####


# Collects single sample files from different locations by creating symlinks
# The output rule only specifies a directory and not the individual symlinks.
# The directory is later symlinked into the working directory of the different
# workflows (see link_input, which is a checkpoint rule, so sample names are
# accessible to snakemake from there on.
# The output directory is named "demux", while samples that are newly demultiplexed
# are placed in other directories (demux_<program>_<method>).
# TODO: demuxing not yet implemented
rule collect_sample_files:
    input:
        sample_tab=run_config + "/samples.tsv"
    output:
        sample_dir=directory("input/{technology}/{layout}/demux/{run}"),
    log:
        "logs/prepare/collect_samples/{technology}/demux/{layout}_{run}.log",
    wildcard_constraints:
        technology = r"\w+",
        layout = r"(single|paired)",
        run = r"\w+",
    group:
        "run"
    script:
        "../scripts/collect_sample_files.py"


# Obtain a list of files from different runs to pool per sample (if pool_raw is true).
# This can only be done within runs with the same demultiplexing method.
# The method collects the file paths to pool into a YAML file
# and creates a samples.tsv file.
rule make_pooling_list:
    input:
        sample_files=lambda wildcards: cfg.get_run_data(
            run=wildcards.run_list + "_pool",
            layout=wildcards.layout
        )["sample_files"],
    output:
        yml="input/sample_config/{technology}/{layout}/{run_list}_pool/samples.yml",
        sample_file="input/sample_config/{technology}/{layout}/{run_list}_pool/samples.tsv",
    log:
        "logs/input/{technology}/{layout}/{run_list}_make_pooling_list.log",
    conda:
        "envs/basic.yaml"
    script:
        "../scripts/make_pooling_list.py"


# Pool samples from different runs (if pool_raw is true).
# This method will just symlink sample files that are only found in a single run,
# meaning that for runs with non-overlapping sample names, no actual pooling of
# sequencing reads is done, however all samples will be treated as coming from
# the same run in the downstream analysis
# (which may be problematic with pipelines such as DADA2)
rule pool_runs_raw:
    input:
        yml=run_config + "_pool/samples.yml",
    output:
        fq=directory(rules.collect_sample_files.output.sample_dir + "_pool"),
    log:
        "logs/input/{technology}/{layout}/demux/{run}_pool_raw.log",
    wildcard_constraints:
        technology = r"\w+",
        layout = r"(single|paired)",
        run = r"\w+",
    threads: workflow.cores,
    conda:
        "envs/basic.yaml"
    script:
        "../scripts/pool_raw.py"


# Symlinks run directories from input to processing/{workflow}/input, selecting
# the ones that were demultiplexed as configured in the workflow.
# *note*: we actually create a "nested" directory, which is used then further.
# The reason is that .snakemake_timestamp in the source directory will interfere
# with the timestamp in the target directory, leading to a lot of weird
# re-running of the workflows
rule link_input:
    input:
        run_dir=rules.collect_sample_files.output.sample_dir,
        # with demultiplexing implemented:
        # run_dir=lambda w: directory(expand("input/{{technology}}/{{layout}}/{demux_method}/{{run}}", demux_method=cfg[w.workflow]["demux_method"])),
    output:
        sample_dir=directory("processing/{workflow}/input/{technology}/{layout}/{run}"),
    wildcard_constraints:
        technology = r"\w+",
        layout = r"(single|paired)",
        run = r"\w+",
    shell:
        """
        mkdir -p {output}
        ln -sr {input.run_dir} {output}/nested
        """


checkpoint final_sample_tab:
    params:
        # we actually have samples in input.sample_dir / subdir
        subdir="nested"
    input:
        tab=run_config + "/samples.tsv",
        sample_dir=rules.link_input.output.sample_dir
    output:
        tab="processing/{workflow}/input/sample_config/{technology}/{layout}/{run}/samples.tsv",
    log:
        "logs/{workflow}/prepare/{technology}_{run}_{layout}/make_final_sample_tab.log",
    script:
        "../scripts/make_new_sample_tab.py"


rule list_samples_yaml:
    input:
        sample_files=lambda wildcards: [
            (run_config + "/samples.yaml").format(**d)
            for d in cfg.get_runs(wildcards.workflow, pooled=False)
        ],
    output:
        yml="results/{workflow}/samples.yaml",
    log:
        "logs/{workflow}/prepare/list_samples_yaml.log",
    script:
        "../scripts/list_samples_yaml.py"


rule collect_unique_files:
    params:
        run_meta=lambda _: list(cfg.get_runs(pooled=False)),
        path_template=lambda _: run_config + "/samples.tsv",
    input:
        sample_files=[
            (run_config + "/samples.tsv").format(**d)
            for d in cfg.get_runs(pooled=False)
        ],
        # fq=lambda wildcards: expand_runs(
        #     path="input/{technology}/{layout}/{demux}/{run}",
        #     **wildcards
        # )
    output:
        read_dir=directory("unique_samples/demux/read_files"),
        tsv="unique_samples/demux/samples.tsv",
    log:
        "logs/unique_samples.log",
    script:
        "../scripts/collect_unique_files.py"



#### Processing ####


# Creates a directory called 'data' inside the workflow results dir 
# if there is only a single run/layout/primer combination in the whole datataset,
# so it is not necessary to dive into all the nested directories
# to get to the results.
# The .outdirs file is actually only used as a flag to ensure that the rule
# is executed. Just deleting the data directory symlink will not lead to this rule
# being re-run, unfortunately.
rule link_data_dir:
    params:
        data_dir="results/{workflow}/data",
    input:
        clust=lambda wildcards: result_paths("/denoised.fasta", workflows=[wildcards.workflow]),
    output:
        touch("results/{workflow}/.outdirs"),
    log:
        "logs/{workflow}/link_data_dir.log",
    script:
        "../scripts/link_data_dir.py"



#### Steps after clustering ####


rule otutab_to_biom:
    input:
        otutab="{prefix}/denoised_otutab.txt.gz",
    output:
        tmp_tab=temp("{prefix}/denoised_otutab_tmp.txt"),
        biom="{prefix}/denoised.biom",
    log:
        "logs/convert_biom/biom_{prefix}.log",
    group:
        "denoise"
    priority: -100
    conda:
        "envs/biom.yaml"
    shell:
        """
        gzip -dc {input.otutab} > {output.tmp_tab}
        biom convert -i {input.otutab} \
          -o {output.biom} \
          --table-type 'OTU table' --to-json &> {log}
        """

rule biom_to_hdf5:
    input:
        biom="{prefix}/denoised.biom",
    output:
        biom_hdf5="{prefix}/denoised.hdf5.biom",
    log:
        "logs/convert_biom/hdf5_{prefix}.log",
    group:
        "denoise"
    priority: -100
    conda:
        "envs/biom.yaml"
    shell:
        # biom convert -i {input.tab} \
        #   -o {output.biom} \
        #   --table-type 'OTU table' --to-json &> {log}
        """
        biom convert -i {input.biom}  \
          -o {output.biom_hdf5} \
          --table-type "OTU table" --to-hdf5 &> {log}
        """


rule combine_sample_reports:
    input:
        reports=lambda wildcards: run_results("/sample_report.tsv", workflows=[wildcards.workflow], allow_missing=True),
    output:
        report="results/{workflow}/sample_report.tsv"
    log:
        "logs/{workflow}/combine_sample_reports.log"
    wildcard_constraints:
        workflow = r"\w+",
        technology = r"\w+",
        layout = r"(single|paired)",
        run = r"\w+",
    script:
        "../scripts/combine_sample_reports.py"

