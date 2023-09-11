from os.path import dirname
import os

from snakemake.workflow import srcdir

from lib import Config

# Set environment variable of workflow root dir to allow post-deploy
# scripts to run other scripts stored in that directory
os.environ['PIPELINE_DIR'] = dirname(dirname(dirname(srcdir('.'))))


#### Configuration ####

cfg = Config(config)

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


def iter_runs(workflows=cfg.workflows, pooled=True):
    for workflow in workflows:
        for run_data in cfg.get_runs(workflow, pooled=pooled):
            wcfg = cfg.workflows[workflow]
            if (run_data["technology"], run_data["layout"]) in cfg.cluster_capabilities[wcfg["pipeline"]]:
                # print(workflow, wcfg["cluster"], run_data)
                yield workflow, wcfg, run_data


def expand_runs(path, workflows=cfg.workflows, pooled=True, **param):
    for workflow, wcfg, run_data in iter_runs(workflows, pooled):
        yield from expand(
            path,
            workflow=workflow,
            **wcfg,
            **run_data,
            **param
        )

def run_results(sub_path="", workflows=cfg.workflows, **param):
    return expand_runs(
        run_result_path + "{sub_path}",
        workflows=workflows,
        sub_path=sub_path,
        **param
    )

# assists in listing results files
# requires sub-path in the results dir, or a function that accepts all workflow settings and returns a sub-path
def result_paths(sub_path="", workflows=cfg.workflows, pooled=True):
    for workflow, wcfg, run_data in iter_runs(workflows, pooled):
        for marker, primer_comb in cfg.primer_combinations.items():
            yield from expand(
                run_result_path + "/{marker}__{primers}{sub_path}",
                workflow=workflow,
                marker=marker,
                primers=primer_comb,
                sub_path=sub_path(marker, primer_comb, wcfg) if callable(sub_path) else sub_path,
                **wcfg,
                **run_data
            )


def tax_dirs(sub_path, ext, workflows=cfg.workflows, pooled=True):
    return result_paths(lambda marker, _, p: expand(
        "/{sub_path}/{name}.{ext}",
        sub_path=sub_path,
        ext=ext,
        name=[
            "{db_name}-{assign[classifier]}-{method_name}".format(**settings)
            for settings in list_taxonomy_runs(p["name"], marker)
        ]
    ), workflows, pooled)


def list_taxonomy_runs(workflow, marker):
    """
    Lists taxonomy settings, removing impossible combinations
    (preformatted/trained databases that cannot be used as input for a given method).
    We don't do the filtering earlier in the Config object because
    cfg.taxonomy_formats is only known after all snakefiles are parsed, while
    the Config object is constructed earlier.
    * training settings are forced to 'standard' (empty) if dealing with a database
       that cannot be imported, but has to be used as-is (pre-formatted/trained)
    * the 'preformatted' option is also set here
    """
    marker_db_cfg = cfg.workflows[workflow]["taxonomy"][marker]
    for d in marker_db_cfg.values():
        db_source_cfg = d["source"]
        classifier = d["assign"]["classifier"]
        expected_format = cfg.taxonomy_formats[classifier]
        if db_source_cfg['format'] in cfg.imported_tax_formats:
            db_source_cfg['preformatted'] = 'regular'
        else:
            db_source_cfg['preformatted'] = 'preformatted'
            d["assign"]["train"] = {"conversion_id": "standard"}
        if db_source_cfg["preformatted"] == "regular" or db_source_cfg["format"] == expected_format:
            yield d


def get_refdb_path(format, ext="", **settings):
    settings = cfg.tax_config(**settings)
    return "refdb/taxonomy/db_{source[preformatted]}_{source[source_id]}/flt_{filter[filter_id]}/{format}/cnv_{assign[train][conversion_id]}{ext}".format(
        format=format, ext=ext, **settings
    )
