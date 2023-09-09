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


#### Configuration ####

cfg = lib.Config(config)

# commonly used
run_config_path = "input/sample_config/{technology}/{layout}/{run}"
run_result_path = "results/{workflow}/workflow_{cluster}/{run}_{layout}"


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
        sorted(set(w.read)) == d["read"], (
        "Sample tab does not align with actual files, "
        "try deleting workdir/<workflow>/input and re-run")
    return expand(full_path if path is None else path, **wildcards, **d)


def expand_workflows(path, workflows=cfg.workflows, **param):
    for workflow in workflows:
        yield from expand(
            path,
            workflow=workflow,
            **cfg.workflows[workflow],
            **param
        )


def expand_runs(path, workflows=cfg.workflows, pooled=True, **param):
    for workflow in workflows:
        for run_data in cfg.get_runs(workflow, pooled=pooled):
            yield from expand(
                path,
                workflow=workflow,
                **cfg.workflows[workflow],
                **run_data,
                **param
            )


def run_results(sub_path="", workflows=cfg.workflows, **param):
    return expand_runs(run_result_path + "{sub_path}", workflows=workflows, sub_path=sub_path, **param)


# assists in listing results files
# requires sub-path in the results dir, or a function that accepts all workflow settings and returns a sub-path
# TODO: redundancy with run_results
def result_paths(sub_path="", workflows=cfg.workflows, pooled=True, allow_missing=True):
    for workflow in workflows:
        p = cfg.workflows[workflow]
        for marker, primer_comb in cfg.primer_combinations.items():
            for r in cfg.get_runs(workflow, pooled=pooled):
                yield from expand(
                    run_result_path + "/{marker}__{primers}{sub_path}",
                    workflow=workflow,
                    cluster=p["cluster"],
                    marker=marker,
                    primers=primer_comb,
                    sub_path=sub_path(marker, primer_comb, p) if callable(sub_path) else sub_path,
                    allow_missing=allow_missing,
                    **r
                )
