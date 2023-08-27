from os.path import *
import shutil
import lib
import os
import csv
import re
from collections import OrderedDict, defaultdict
from glob import glob
import copy


localrules:
    get_run_config,
    make_pooling_list,
    dump_config,
    make_outdirs,
    link_data_dir,
    combine_sample_reports,

ruleorder:
    # if pool_raw is true, an artificial run called 'pool' will be created,
    # resulting in an ambiguous situation
    pool_runs_raw > collect_samples


#### Start ####

cfg = lib.Config(config)


#### Helper functions ####


def cfg_fail(msg):
    # cfg_dir = os.path.join("input", "sample_config")
    # if os.path.exists(cfg_dir):
    #     shutil.rmtree(cfg_dir)
    raise Exception(msg + " Try re-running.")


def get_runs(technology="*", layout="*", run="*", pool=False):
    # print("get runs", technology, layout, run, pool)
    # ensure that 'get_run_config' completed
    cfg_out = checkpoints.get_run_config.get().output[0]
    # then, expand layout subdirectories, either by inserting the glob star (*)
    # or some custom supplied wildcards
    exp = expand("{cfg_out}/{technology}/{layout}",
                 cfg_out=cfg_out,
                 technology=technology,
                 layout=layout)
    for pattern in exp:
        g = glob(pattern)
        if len(g) == 0:
            cfg_fail("Configuration incomplete.")
        for layout_dir in g:
            p = layout_dir.split('/')
            par = {"technology": p[-2], "layout": p[-1]}
            if run.endswith("_pool"):
                # the pooled sample file may not exist yet, so we can't use glob
                assert pool, f"get_runs: pool must be True with run {run}"
                par["run"] = run
                yield os.path.join(layout_dir, run), copy.copy(par)
            else:
                # get runs in similar way: first use expand(), then glob()
                runs = []
                for pattern in expand("{layout_dir}/{run}", layout_dir=layout_dir, run=run):
                    for run_dir in glob(pattern):
                        par["run"] = run_dir.split('/')[-1]
                        # "pool" is special and should be excluded
                        if not par["run"].endswith("_pool"):
                            runs.append((run_dir, copy.copy(par)))
                if len(runs) == 0:
                    cfg_fail(f"Config directory empty: {layout_dir}.")
                if pool and len(runs) > 1:
                    # The pooled run name is formed from a sorted list of
                    # run names. This is the only place where the this name
                    # is defined. At other places in the code, we only need
                    # to know (and enforce) that runs ending with "_pool" are
                    # in fact run pools.
                    run_names = sorted(os.path.basename(r) for r, _ in runs)
                    par["run"] = "_".join(run_names) + "_pool"
                    yield os.path.join(layout_dir, par["run"]), par
                else:
                    yield from iter(runs)



def expand_runs(path, technology="*", layout="*", run="*", pool=False, **params):
    """
    Expand-like function, additionally providing run=<run list>
    """
    out = []
    for _, p in get_runs(technology=technology, layout=layout, run=run, pool=pool):
        params.update(p)
        # expand supplied wildcards other than the ones known by get_runs() (technology, layout, run)
        for exp_path in expand(path, **params):
            out.append(exp_path)
    assert len(set(out)) == len(out), "Duplicate output files, some wildcards need a value"
    return out


sample_re = re.compile("(.+?)_R1\.fastq\.gz")


def expand_samples(path=None, technology="*", layout="*", run="*", pool=False, **params):
    """
    Expand-like function, additionally providing run=<run list> and
    sample=<sample list> based on dynamically executed checkpoint rules.
    The 'technology' and 'layout' wildcards can be specified for obtaining samples
    only for specific parameter combinations. If not, all combinations will be
    returned.
    """
    outfiles = []
    for run_dir, p in get_runs(technology=technology, layout=layout, run=run, pool=pool):
        # ensure that 'collect_samples', ('combine_runs') and 'link_input' have
        # been executed
        params.update(p)
        out2 = checkpoints.link_input.get(**params).output[0]
        # there are two ways of obtaining the samples:
        # (1) from the sample file
        sample_file = os.path.join(run_dir, "samples.tsv")
        if not os.path.exists(sample_file):
            sample_file = os.path.join(run_dir, "indexes.tsv")
            if not os.path.exists(sample_file):
                cfg_fail(f"Neither samples.tsv nor indexes.tsv present in {run_dir}.")
        with open(sample_file) as f:
            rdr = csv.reader(f, delimiter='\t')
            next(rdr)  # skip header
            samples = sorted(row[0] for row in rdr if row)
        # (2) by using glob
        samples2 = sorted(sample_re.fullmatch(s.split('/')[-1]).group(1)
                         for s in glob(out2 + "/*_R1.fastq.gz"))
        # the results should be the same
        assert samples == samples2, "Sample file mismatch: delete the 'input' and/or 'processing' directory and retry"
        reads = [1, 2] if p["layout"] == "paired" else [1]
        outfiles += expand(path, sample=samples, read=reads, **params)
    return outfiles


#### Prepare / configure ####


# Obtains sample file lists ("manifest files") for every run and run layout
# given base directories, file patterns or already assembled manifest files
checkpoint get_run_config:
    params:
        cfg = config["input"]
    output:
        dir=directory("input/sample_config")
    log:
        "logs/prepare/get_run_config.log",
    script:
        "../scripts/get_run_config.py"


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
#         indexes="input/sample_config/{technology}/{layout}/{run}/runfiles.tsv"
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
#         indexes="input/sample_config/{technology}/{layout}/{run}/indexes.tsv"
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
rule collect_samples:
    input:
        samples="input/sample_config/{technology}/{layout}/{run}/samples.tsv"
    output:
        fq=directory("input/{technology}/{layout}/demux/{run}"),
    log:
        "logs/prepare/collect_samples/{technology}/{layout}/{run}.log",
    group:
        "run"
    script:
        "../scripts/collect_sample_files.py"


# Obtain a list of files from different runs to pool per sample (if pool_raw is true).
# This can only be done within runs with the same demultiplexing method.
# The method collects the file paths to pool into a YAML file
# and creates a samples.tsv file.
# TODO: For runs not yet demultiplexed, sample names from indexes.tsv are used to construct file names where expected after demultiplexed
rule make_pooling_list:
    input:
        # run_config=lambda wildcards: expand_runs(
        #     "input/sample_config/{technology}/{layout}/{run}",
        #     **wildcards
        # ),
        run_config=lambda wildcards: expand_runs(
            "input/sample_config/{technology}/{layout}/{run}",
            **wildcards
        ),
    output:
        directory("input/sample_config/{technology}/{layout}/{run_list}_pool"),
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
        yml="input/sample_config/{technology}/{layout}/{run_list}_pool/samples.yml",
    output:
        fq=directory("input/{technology}/{layout}/{demux_method}/{run_list}_pool"),
    log:
        "logs/input/{technology}/{layout}/{demux_method}/{run_list}_pool_raw.log",
    threads: workflow.cores,
    conda:
        "envs/basic.yaml"
    script:
        "../scripts/pool_raw.py"


# Symlinks run directories from input to processing/{workflow}/input, selecting
# the ones that were demultiplexed as configured in the workflow.
# This is a checkpoint: sample names are only known to snakemake after
#  execution of this rule.
checkpoint link_input:
    input:
        "input/{technology}/{layout}/demux/{run}",
        # with demultiplexing implemented:
        # fq=lambda w: directory(expand("input/{{technology}}/{{layout}}/{demux_method}/{{run}}", demux_method=cfg[w.workflow]["demux_method"])),
    output:
        directory("processing/{workflow}/input/{technology}/{layout}/{run}"),
    wildcard_constraints:
        technology = r"\w+",
        layout = r"(single|paired)",
        run = r"\w+",
    shell:
        """
        ln -sr {input} {output}
        """


rule list_samples_yaml:
    input:
        sample_files=lambda wildcards: expand_runs(
            "input/sample_config/{technology}/{layout}/{run}/samples.yaml",
            **wildcards
        ),
    output:
        yml="results/{workflow}/samples.yaml",
    log:
        "logs/{workflow}/prepare/list_samples_yaml.log",
    script:
        "../scripts/list_samples_yaml.py"


rule collect_unique_files:
    input:
        run_cfg_dirs=lambda wildcards: expand_runs(
            path="input/sample_config/{technology}/{layout}/{run}",
            **wildcards
        ),
        # fq=lambda wildcards: expand_runs(
        #     path="input/{technology}/{layout}/demux/{run}",
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

# Creates the output directories (needed for link_data_dir to work)
rule make_outdirs:
    priority: 100
    output:
        touch(directory("results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{fprimer}...{rprimer}"))
    wildcard_constraints:
        layout = r"(single|paired)",
        run = r"\w+",
        cluster = r"\w+",
        marker = r"\w+",
        fprimer = r"[^/]+",
        rprimer = r"[^/]+",


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
        outdirs=lambda wildcards: expand_runs(
            "results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            cluster=cfg[wildcards.workflow]["cluster"],
            primers=cfg.primer_combinations_flat,
            **wildcards
        )
    output:
        list="results/{workflow}/.outdirs",
    log:
        "logs/{workflow}/link_data_dir.log",
    script:
        "../scripts/link_data_dir.py"



#### Steps after clustering ####


rule tsv_to_biom:
    input:
        tab="results/{workflow}/workflow_{cluster}/{primers}/{layout}_{run}/denoised_otutab.txt.gz",
    output:
        biom="results/{workflow}/workflow_{cluster}/{primers}/{layout}_{run}/denoised.biom",
        biom_hdf5="results/{workflow}/workflow_{cluster}/{primers}/{layout}_{run}/denoised.hdf5.biom",
    log:
        "logs/{workflow}/workflow_{cluster}/{primers}/{layout}_{run}/tsv_to_biom.log",
    group:
        "denoise"
    conda:
        "envs/biom.yaml"
    shell:
        """
        biom convert -i {input.tab} \
          -o {output.biom} \
          --table-type 'OTU table' --to-json &> {log}
        biom convert -i {output.biom}  \
          -o {output.biom_hdf5} \
          --table-type "OTU table" --to-hdf5 &> {log}
        """


rule combine_sample_reports:
    input:
        reports=lambda wildcards: chain(*[expand_runs(
            "results/{workflow}/workflow_{cluster}/{run}_{layout}/sample_report.tsv",
            cluster=cfg[wildcards.workflow]["cluster"],
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards
        )])
    output:
        report="results/{workflow}/sample_report.tsv"
    log:
        "logs/{workflow}/combine_sample_reports.log"
    script:
        "../scripts/combine_sample_reports.py"
