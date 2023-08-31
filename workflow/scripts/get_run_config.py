import glob
import os
import shutil

from utils import file_logging
from utils.make_manifest import make_manifest
from utils.sample_list import SampleList


def get_config_file(outdir, run_name, sample_file=None, **param):
    reserved_chars = " -."
    runs = []
    if sample_file is None:
        # search for read files with given parameters
        _runs = make_manifest(
            out_prefix=outdir,
            path_template="{out_prefix}/{layout}/{run}/samples.tsv",
            header_single=SampleList.default_header["single"],
            header_paired=SampleList.default_header["paired"],
            default_run=run_name, 
            relative_paths=True,
            reserved_chars=reserved_chars,
            **param
        )
        runs.extend(r for r, _ in _runs)
    else:
        # sample list supplied -> read and normalize sample names
        sample_list = SampleList(sample_file, reserved_chars=reserved_chars)        
        # write to new location
        _outdir = os.path.join(outdir, sample_list.layout, run_name)
        outfile = os.path.join(_outdir, "samples.tsv")
        if not os.path.exists(_outdir):
            os.makedirs(_outdir)
        with open(outfile, "w") as o:
            sample_list.write(o)
        runs.append(run_name)
    # read sample files again in order to write it to YAML
    for f in glob.glob(f"{outdir}/*/*/samples.tsv"):
        sample_list = SampleList(f)
        yml_out = os.path.join(os.path.dirname(f), "samples.yaml")
        with open(yml_out, "w") as o:
            sample_list.write_yaml(o)
    return runs


def get_run_config(outdir, cfg):
    # important: remove the whole output directory before adding anything, since
    # the checkpoint rules depend on the contents of this dir
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    is_entry = lambda cfg: isinstance(cfg, dict) and (
        'sample_pattern' in cfg and ("directories" in cfg or "patterns" in cfg) \
        or "sample_file" in cfg
    )
    assert isinstance(cfg, dict)
    runs = []
    for technology, tcfg in cfg.items():
        _outdir = os.path.join(outdir, technology)
        if is_entry(tcfg):
            runs += get_config_file(_outdir, "run1", **tcfg)
        else:
            for run, run_cfg in tcfg.items():
                assert is_entry(run_cfg), (
                    "Invalid 'input' configuration: {{ {}: {} }}. Either specify "
                    "the input files for a single run directly below 'input' "
                    "(sample_pattern, directories/patterns, ...). "
                    "Or add a configuration for multiple runs, each with "
                    "these separate keywords. The third option is to use the "
                    "{{run}} wildcard in directories/paths to obtain files from "
                    "multiple runs.").format(run, str(run_cfg))
                rw = "{run}"
                assert not rw in run and not any(rw in d for d in run_cfg.get("directories", []) + run_cfg.get("patterns", [])), (
                    "No {run} wildcards allowed in run names or directories/patterns"
                    " of named runs")
                runs += get_config_file(_outdir, run, **run_cfg)
    # check for duplicate run names
    assert len(set(runs)) == len(runs), (
        "Duplicate run name(s) detected. Runs must be unique even "
        "across sequencing technologies")
    for r in runs:
        assert not r.endswith("_pool"), \
            f"Run '{r}' should not end with '_pool' (reserved for run pools)"


with file_logging(snakemake.log[0]):
    get_run_config(snakemake.output.dir, snakemake.params.cfg)
