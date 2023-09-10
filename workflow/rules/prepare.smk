
localrules:
    collect_sample_lists,
    make_pooling_list,
    dump_config,
    link_data_dir,
    collect_unique_files


ruleorder:
    # If pool_raw is true, an artificial run called 'pool' will be created,
    # resulting in an ambiguous situation, where run pooling needs to have priority
    pool_runs_raw > collect_sample_files


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
        path_template = lambda _: run_config_path + "/samples.{ext}",
    output:
        sample_files=sorted(chain(*[
            expand(
                run_config_path + "/samples.{ext}",
                ext=["tsv", "yaml"],
                **d
            )
            for d in cfg.get_runs(pooled=False)
        ])),
    log:
        "logs/input/collect_sample_lists.log",
    script:
        "../scripts/collect_sample_lists.py"


rule dump_config:
    params:
        config=lambda wildcards: cfg.workflows[wildcards.workflow],
    output:
        "results/{workflow}/config.yaml",
    log:
        "logs/{workflow}/dump_config.log",
    script:
        "../scripts/dump_config.py"


# checks and reverse complements primer sequences and calculates the consensus
# of primer mixtures
rule prepare_primers:
    params:
        primers=cfg.primers,
    output:
        yaml="workdir/primers/primers.yaml",
    log:
        "logs/prepare_primers.log",
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
#         indexes=run_config_path + "/indexes.tsv"
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
#         indexes=run_config_path + "/indexes.tsv"
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
        sample_tab=run_config_path + "/samples.tsv"
    output:
        sample_dir=directory("input/{technology}/{layout}/demux/{run}"),
    log:
        "logs/input/{technology}/{layout}/{run}/collect_sample_files.log",
    wildcard_constraints:
        technology = r"\w+",
        layout = r"(single(\.rev)?|paired)",
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
        "logs/input/{technology}/{layout}/{run_list}_pool/make_pooling_list.log",
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
        yml=run_config_path + "_pool/samples.yml",
    output:
        fq=directory(rules.collect_sample_files.output.sample_dir + "_pool"),
    log:
        "logs/input/{technology}/{layout}/{run}_pool/pool_raw.log",
    wildcard_constraints:
        technology = r"\w+",
        layout = r"(single(\.rev)?|paired)",
        run = r"\w+",
    threads: workflow.cores,
    conda:
        "envs/basic.yaml"
    script:
        "../scripts/pool_raw.py"


# Symlinks run directories from input to workdir/{workflow}/input, selecting
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
        sample_dir=directory("workdir/{workflow}/input/{technology}/{layout}/{run}"),
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
        tab=run_config_path + "/samples.tsv",
        sample_dir=rules.link_input.output.sample_dir
    output:
        tab="workdir/{workflow}/input/sample_config/{technology}/{layout}/{run}/samples.tsv",
    log:
        "logs/{workflow}/{run}_{layout}/{technology}_make_final_sample_tab.log",
    script:
        "../scripts/make_new_sample_tab.py"


rule list_samples_yaml:
    input:
        sample_files=lambda wildcards: [
            (run_config_path + "/samples.yaml").format(**d)
            for d in cfg.get_runs(wildcards.workflow, pooled=False)
        ],
    output:
        yml="results/{workflow}/samples.yaml",
    log:
        "logs/{workflow}/list_samples_yaml.log",
    script:
        "../scripts/list_samples_yaml.py"


rule collect_unique_files:
    params:
        run_meta=lambda _: list(cfg.get_runs(pooled=False)),
        path_template=lambda _: run_config_path + "/samples.tsv",
    input:
        sample_files=[
            (run_config_path + "/samples.tsv").format(**d)
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



#### Workdir ####


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
        clust=lambda wildcards: result_paths("/clusters.fasta", workflows=[wildcards.workflow]),
    output:
        out_list="results/{workflow}/.outdirs",
    log:
        "logs/{workflow}/link_data_dir.log",
    script:
        "../scripts/link_data_dir.py"

