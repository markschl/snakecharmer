import os
from os.path import exists, dirname, join, relpath

from utils import file_logging


def link_samples(results_dirs, workflow_cfg):
    for sym_dir in results_dirs:
        p = workflow_cfg[os.path.basename(dirname(sym_dir))]
        res_dir = join(
            dirname(sym_dir),
            "pipeline_" + p["cluster"],
            p["single_primercomb"],
            p["single_strategy"],
        )
        if not exists(res_dir):
            os.makedirs(res_dir)
        if exists(sym_dir):
            os.remove(sym_dir)
        os.symlink(relpath(res_dir, os.path.dirname(sym_dir)), sym_dir)


with file_logging(snakemake.log[0]):
    link_samples(results_dirs=snakemake.output, workflow_cfg=snakemake.params.workflow_cfg)
