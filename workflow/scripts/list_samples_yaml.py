from collections import defaultdict
import os

import yaml

from utils import file_logging


def list_samples(sample_files, outfile):
    out = defaultdict(dict)
    for sample_file in sample_files:
        p = sample_file.split(os.sep)
        layout = p[-3]
        run = p[-2]
        with open(sample_file) as f:
            out[run][layout] = yaml.safe_load(f)
    with open(outfile, "w") as o:
        yaml.safe_dump(dict(out), o)


with file_logging(snakemake.log[0]):
    list_samples(sample_files=snakemake.input.sample_files,
                    outfile=snakemake.output.yml)
