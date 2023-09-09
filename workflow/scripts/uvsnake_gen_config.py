from os.path import relpath, dirname

import yaml

from utils import file_logging


snakefile_content = """
configfile: "{configfile}"

module uvsnake:
    snakefile: "{snakefile}"
    config: config

use rule * from uvsnake as uvsnake_*

# special target rules for preprocessing
# before applying either UNOISE3 or UPARSE
rule uvsnake_prepare:
    input:
        trim_dir=expand(
            "workdir/prepare_paired/2_trim/{{sample}}/{{sample}}_{{dir_}}.log",
            sample=uvsnake.config["_sample_names"],
            dir_=["fwd", "rev"]
        ),
        report="results/sample_report.tsv",
"""


def write_config(sample_tab, config_out, snakefile, snakefile_out, primer_config, usearch_config):
    # generate the configuration
    out = {}
    workdir = dirname(snakefile_out)
    out["input"] = {"sample_file": relpath(sample_tab, workdir)}
    # prepare primers:
    # uvsnake has almost the same configuration,  but does not have the
    # marker concept of, so we merge primers from all markers.
    # Primers have been ensured to be unique across markers, so
    # there will not be any name clashes
    out["primers"] = {
        "forward": [],
        "reverse": [],
        "trim_settings": primer_config.pop("trim_settings")
    }
    for cfg in primer_config.values():
        for _dir, primers in cfg.items():
            out["primers"][_dir] += primers

    config_keys = ["defaults", "merge", "filter", "unoise3", "uparse", "otutab"]
    for k in config_keys:
        out[k] = usearch_config[k]
    out["merge"].pop("expected_length", None)

    with open(config_out, "w") as f:
        yaml.safe_dump(out, f, sort_keys=False)

    # generate the Snakefile
    with open(snakefile_out, "w") as o:
        o.write(snakefile_content.format(
            configfile=relpath(config_out, workdir),
            snakefile=snakefile,
        ))
    

with file_logging(snakemake.log[0]):
    write_config(
        sample_tab=snakemake.input.sample_tab,
        config_out=snakemake.output.config,
        snakefile_out=snakemake.output.snakefile,
        **snakemake.params
    )
